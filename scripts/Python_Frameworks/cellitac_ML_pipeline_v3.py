#!/usr/bin/env python3

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier

from sklearn.model_selection import (train_test_split, StratifiedKFold,
                                     cross_val_predict, learning_curve)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import (classification_report, confusion_matrix,
                             accuracy_score, precision_score, recall_score,
                             f1_score, balanced_accuracy_score, roc_auc_score)

from imblearn.pipeline import Pipeline as ImbPipeline
from imblearn.over_sampling import SMOTE

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid")
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# ----------------------------------------------------------------------------- 
# CONFIGURATION
# -----------------------------------------------------------------------------
MIN_CELLS        = 50     # classes with fewer cells cannot be cross-validated honestly
MAX_FEATURES     = 150    # cap on motifs entering the classifier (selected in-fold)
USE_SMOTE        = False  # False -> class_weight="balanced" (no synthetic cells)
CORR_CUTOFF      = 0.95   # unsupervised redundancy filter (uses no labels)

# --- statistical criterion for TF <-> cell-type association (# NEW) ----------
FDR_ALPHA        = 0.05   # Benjamini-Hochberg level across all TF x class tests
MIN_COHENS_D     = 0.50   # minimum standardised effect size (medium effect)
MIN_RB_AUC       = 0.60   # minimum rank-biserial AUC (0.5 = no separation)

# --- permutation importance (only needed for models without native importance)
PERM_N_REPEATS   = 5
PERM_SUBSAMPLE   = 1000   # test cells used per permutation pass (speed)

# --- network drawing ---------------------------------------------------------
NET_TOP_TFS      = 20     # candidate TFs per model, ranked by that model
SHOW_TF_TF_EDGES = False  # off by default: the requested graph is cell <-> TF

# --- learning curve ----------------------------------------------------------
LC_TRAIN_SIZES   = np.linspace(0.1, 1.0, 6)
LC_CV_FOLDS      = 3

# RNA markers that must NEVER appear as ATAC features (leakage tripwire)
RNA_MARKERS = {"GNLY", "NKG7", "PRF1", "SLC8A1", "PLXDC2", "CD3D", "CD3E",
               "CD8A", "CD4", "MS4A1", "LYZ", "CD14", "FCGR3A", "PPBP",
               "BLNK", "BCL11A", "MARCH1", "CD19", "NCAM1"}

PALETTE = "tab10"


# =============================================================================
class CellitacPipeline:

    def __init__(self, data_dir="cellitac_output", 
                 out_dir="output/output_Second_Framework_after_dropping"):
        self.data_dir = Path(data_dir)
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.models, self.fitted, self.results = {}, {}, {}
        self.label_encoder = LabelEncoder()
        print("=" * 78)
        print("cellitac pipeline v3  |  leakage-safe, TF-centric")
        print(f"  data dir : {self.data_dir.resolve()}")
        print(f"  out dir  : {self.out_dir.resolve()}")
        print("=" * 78)

    def _save(self, fig, name):
        path = self.out_dir / name
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"      [fig] {path.name}")

    def _save_table(self, df, name, title=None, float_fmt="{:.4f}"):
        csv_path = self.out_dir / f"{name}.csv"
        df.to_csv(csv_path, index=False)
        disp = df.copy()
        for c in disp.columns:
            if pd.api.types.is_float_dtype(disp[c]):
                disp[c] = disp[c].map(lambda v: float_fmt.format(v) if pd.notna(v) else "")
        n_rows, n_cols = disp.shape
        fig, ax = plt.subplots(figsize=(min(2.0 * n_cols, 20), 0.55 * n_rows + 1.4))
        ax.axis("off")
        tbl = ax.table(cellText=disp.values, colLabels=disp.columns,
                       cellLoc="center", loc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(10)
        tbl.scale(1, 1.5)
        for j in range(n_cols):                       # header styling
            cell = tbl[0, j]
            cell.set_facecolor("#4059AD")
            cell.set_text_props(color="white", weight="bold")
        for i in range(1, n_rows + 1):                # zebra striping
            for j in range(n_cols):
                if i % 2 == 0:
                    tbl[i, j].set_facecolor("#EEF1F8")
        if title:
            ax.set_title(title, fontsize=13, weight="bold", pad=18)
        self._save(fig, f"{name}.png")
        print(f"      [tbl] {csv_path.name}")
        return csv_path

    @staticmethod
    def _resolve_overlaps(pos, min_dist=0.16, n_iter=250):

        keys = list(pos)
        P = np.array([pos[k] for k in keys], dtype=float)
        for _ in range(n_iter):
            moved = False
            for i in range(len(P)):
                for j in range(i + 1, len(P)):
                    d = P[i] - P[j]
                    dist = float(np.hypot(*d))
                    if dist < min_dist:
                        if dist < 1e-9:                      # exact overlap
                            d = np.random.default_rng(i * 97 + j).normal(size=2)
                            dist = float(np.hypot(*d))
                        shift = (min_dist - dist) / 2.0 * d / dist
                        P[i] += shift
                        P[j] -= shift
                        moved = True
            if not moved:
                break
        return {k: P[i] for i, k in enumerate(keys)}

    @staticmethod
    def _spread_labels(xy, min_gap):

        pts = np.array(xy, dtype=float)
        order = np.argsort(pts[:, 1])
        for a, b in zip(order[:-1], order[1:]):
            gap = pts[b, 1] - pts[a, 1]
            if gap < min_gap and abs(pts[b, 0] - pts[a, 0]) < min_gap * 2.5:
                pts[b, 1] = pts[a, 1] + min_gap
        return pts

    @staticmethod
    def _bh_fdr(pvals):

        p = np.asarray(pvals, dtype=float)
        n = p.size
        order = np.argsort(p)
        ranked = p[order] * n / (np.arange(n) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity
        out = np.empty(n, dtype=float)
        out[order] = np.clip(ranked, 0, 1)
        return out

    # =========================================================== 1. LOAD DATA
    def load_data(self):
        print("\n[1] Loading chromVAR TF activity (ATAC only) and labels ...")
        X = pd.read_csv(self.data_dir / "cellitac_TF_activity.csv", index_col=0)

        leaked = RNA_MARKERS.intersection(set(map(str, X.columns)))
        if leaked:
            raise ValueError(
                f"RNA-marker columns found in the feature matrix: {sorted(leaked)}. "
                "Expression has leaked into X. Aborting.")

        labels = pd.read_csv(self.data_dir / "cell_labels.csv").set_index("cell_id")
        common = X.index.intersection(labels.index)
        X = X.loc[common]
        y = labels.loc[common, "cell_type"].dropna()
        X = X.loc[y.index]
        print(f"    loaded  : {X.shape[0]} cells x {X.shape[1]} TF motifs")

        # counts BEFORE rare-class removal, kept for the pie chart (# NEW)
        self.counts_before = y.value_counts()

        counts = y.value_counts()
        keep = counts[counts >= MIN_CELLS].index
        dropped = counts[counts < MIN_CELLS]
        if len(dropped):
            print(f"    dropped : classes with < {MIN_CELLS} cells -> {dict(dropped)}")
        mask = y.isin(keep)
        self.X, self.y = X.loc[mask], y.loc[mask]
        self.counts_after = self.y.value_counts()

        self.y_encoded = self.label_encoder.fit_transform(self.y)
        # cell-type names are READ FROM THE FILE, never typed in
        self.class_names = list(self.label_encoder.classes_)
        print(f"    modelled: {self.X.shape[0]} cells | {len(self.class_names)} classes")
        print(f"    classes : {self.class_names}")

        # optional motif id -> TF symbol map, used only for reporting
        self.motif_map = None
        mm = self.data_dir / "motif_to_TF_map.csv"
        if mm.exists():
            self.motif_map = pd.read_csv(mm)
            print(f"    motif map loaded ({len(self.motif_map)} entries)")
        return True

    # ================================================ 2. CLASS COMPOSITION
    def analyze_class_imbalance(self):
        print("\n[2] Class composition ...")
        c = self.counts_after
        self.imbalance_df = pd.DataFrame({
            "Cell_Type": c.index,
            "Count": c.values,
            "Proportion": (c / c.sum()).values,
            "Imbalance_Ratio": (c.max() / c).values,
        })
        print(self.imbalance_df.round(3).to_string(index=False))
        self._save_table(self.imbalance_df, "table01_class_composition",
                         "Class composition of the modelled dataset")

        # ---- fig02: pies before / after rare-class removal (# NEW) ----------
        fig, axes = plt.subplots(1, 2, figsize=(15, 7))
        for ax, counts, sub in [
            (axes[0], self.counts_before,
             f"A. After joint QC (n = {int(self.counts_before.sum()):,} cells, "
             f"{len(self.counts_before)} classes)"),
            (axes[1], self.counts_after,
             f"B. Modelled dataset, classes < {MIN_CELLS} cells removed "
             f"(n = {int(self.counts_after.sum()):,}, {len(self.counts_after)} classes)")]:
            colors = sns.color_palette(PALETTE, len(counts))
            wedges, _, autotexts = ax.pie(
                counts.values, labels=None, autopct=lambda p: f"{p:.1f}%",
                startangle=90, colors=colors, pctdistance=0.78,
                wedgeprops=dict(width=0.55, edgecolor="white", linewidth=1.5))
            for t in autotexts:
                t.set_fontsize(9)
            ax.set_title(sub, fontsize=11, weight="bold")
            ax.legend(wedges, [f"{k}  (n={v:,})" for k, v in counts.items()],
                      loc="center left", bbox_to_anchor=(0.98, 0.5), fontsize=9,
                      frameon=False)
        fig.suptitle("Cell-type composition before and after rare-class removal",
                     fontsize=14, weight="bold")
        fig.tight_layout()
        self._save(fig, "fig02_class_composition_pies.png")

        # ---- fig03: what balancing actually does (# NEW) --------------------
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        order = list(self.counts_after.index)
        colors = sns.color_palette(PALETTE, len(order))

        axes[0].bar(range(len(order)), [self.counts_after[o] for o in order], color=colors)
        axes[0].set_xticks(range(len(order)))
        axes[0].set_xticklabels(order, rotation=35, ha="right", fontsize=9)
        axes[0].set_ylabel("Cells")
        axes[0].set_title(f"A. Observed counts (imbalance ratio "
                          f"{self.imbalance_df['Imbalance_Ratio'].max():.2f}:1)",
                          fontsize=11, weight="bold")
        for i, o in enumerate(order):
            axes[0].text(i, self.counts_after[o], f"{self.counts_after[o]:,}",
                         ha="center", va="bottom", fontsize=8)

        if USE_SMOTE:
            sm = SMOTE(k_neighbors=3, random_state=RANDOM_STATE)
            _, y_res = sm.fit_resample(self.X.values, self.y_encoded)
            res = pd.Series(self.label_encoder.inverse_transform(y_res)).value_counts()
            vals = [res.get(o, 0) for o in order]
            ttl = "B. After SMOTE oversampling (synthetic cells added in-fold)"
            ylab = "Cells after resampling"
        else:
            w = compute_class_weight("balanced", classes=np.unique(self.y_encoded),
                                     y=self.y_encoded)
            wmap = {self.label_encoder.classes_[i]: w[i] for i in range(len(w))}
            vals = [self.counts_after[o] * wmap[o] for o in order]
            ttl = ("B. Effective contribution after class_weight='balanced'\n"
                   "(count x class weight — every class contributes equally)")
            ylab = "Effective weighted cells"
        axes[1].bar(range(len(order)), vals, color=colors)
        axes[1].set_xticks(range(len(order)))
        axes[1].set_xticklabels(order, rotation=35, ha="right", fontsize=9)
        axes[1].set_ylabel(ylab)
        axes[1].set_title(ttl, fontsize=11, weight="bold")
        for i, v in enumerate(vals):
            axes[1].text(i, v, f"{v:,.0f}", ha="center", va="bottom", fontsize=8)
        fig.suptitle("Class-imbalance handling", fontsize=14, weight="bold")
        fig.tight_layout()
        self._save(fig, "fig03_class_imbalance_handling.png")
        return self.imbalance_df

    # ================================== 3. UNSUPERVISED CLEANING (no labels)
    def clean_features(self):
        print("\n[3] Unsupervised feature cleaning (no labels used) ...")
        n0 = self.X.shape[1]
        zero_var = self.X.columns[self.X.var() == 0]
        if len(zero_var):
            self.X = self.X.drop(columns=zero_var)
            print(f"    removed {len(zero_var)} zero-variance motifs")
        corr = self.X.corr().abs()
        upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
        high = [c for c in upper.columns if (upper[c] > CORR_CUTOFF).any()]
        if len(high):
            self.X = self.X.drop(columns=high)
            print(f"    removed {len(high)} redundant motifs (|r| > {CORR_CUTOFF})")
        print(f"    motifs  : {n0} -> {self.X.shape[1]}")
        return self.X

    # ============================== 4. UMAP OF THE INPUT SPACE (before model)
    def plot_input_embedding(self):
        """# NEW — UMAP of the TF-activity space the classifier is handed,
        coloured by the labels read from cell_labels.csv."""
        print("\n[4] UMAP of the TF-activity space (before training) ...")
        Z = StandardScaler().fit_transform(self.X.values)
        try:
            import umap
            emb = umap.UMAP(n_neighbors=30, min_dist=0.3, n_components=2,
                            metric="euclidean", random_state=RANDOM_STATE
                            ).fit_transform(Z)
            method = "UMAP"
        except Exception as e:                      # graceful, never fatal
            print(f"    [warn] umap-learn unavailable ({e}); using PCA + t-SNE")
            emb = TSNE(n_components=2, perplexity=30, init="pca",
                       random_state=RANDOM_STATE).fit_transform(
                           PCA(n_components=50, random_state=RANDOM_STATE).fit_transform(Z))
            method = "t-SNE"
        self.input_embedding, self.input_embedding_method = emb, method

        fig, ax = plt.subplots(figsize=(11, 9))
        colors = sns.color_palette(PALETTE, len(self.class_names))
        anchors = []
        for i, ct in enumerate(self.class_names):
            m = (self.y.values == ct)
            ax.scatter(emb[m, 0], emb[m, 1], s=5, alpha=0.65, color=colors[i],
                       label=f"{ct} (n={int(m.sum()):,})", linewidths=0)
            anchors.append([np.median(emb[m, 0]), np.median(emb[m, 1])])
        anchors = self._spread_labels(anchors, min_gap=0.06 * np.ptp(emb[:, 1]))
        for i, ct in enumerate(self.class_names):
            ax.text(anchors[i, 0], anchors[i, 1], ct, fontsize=10,
                    weight="bold", ha="center", va="center", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec=colors[i], lw=1.5, alpha=0.9))
        ax.set_xlabel(f"{method} 1"); ax.set_ylabel(f"{method} 2")
        ax.set_title(f"{method} of chromVAR TF motif-activity space\n"
                     f"{self.X.shape[0]:,} cells x {self.X.shape[1]} TF motifs — "
                     f"coloured by reference-derived cell type",
                     fontsize=13, weight="bold")
        ax.legend(markerscale=3, fontsize=9, loc="center left",
                  bbox_to_anchor=(1.01, 0.5), frameon=False)
        fig.tight_layout()
        self._save(fig, "fig01_umap_TFactivity_before_training.png")
        return emb

    # ================================================== 5. SPLIT + PIPELINES
    def prepare_train_test_split(self):
        print("\n[5] Stratified split (before any supervised step) ...")
        (self.X_train, self.X_test,
         self.y_train, self.y_test) = train_test_split(
            self.X, self.y_encoded, test_size=0.2,
            random_state=RANDOM_STATE, stratify=self.y_encoded)
        print(f"    train {self.X_train.shape[0]:,} | test {self.X_test.shape[0]:,}")
        return True

    def _make_pipeline(self, clf):
        steps = [("scale", StandardScaler()),
                 ("select", SelectKBest(f_classif,
                                        k=min(MAX_FEATURES, self.X_train.shape[1])))]
        if USE_SMOTE:
            steps.append(("smote", SMOTE(k_neighbors=3, random_state=RANDOM_STATE)))
        steps.append(("clf", clf))
        return ImbPipeline(steps)

    def define_models(self):
        print("\n[6] Defining models ...")
        cw = None if USE_SMOTE else "balanced"
        base = {
            "LogReg": LogisticRegression(penalty="l2", C=0.1, class_weight=cw,
                                         max_iter=3000, random_state=RANDOM_STATE),
            "Random Forest": RandomForestClassifier(
                n_estimators=300, max_depth=6, min_samples_leaf=15,
                max_features="sqrt", class_weight=cw,
                random_state=RANDOM_STATE, n_jobs=-1),
            "XGBoost": XGBClassifier(
                n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8,
                colsample_bytree=0.8, reg_alpha=1.0, reg_lambda=5.0,
                random_state=RANDOM_STATE, eval_metric="mlogloss", verbosity=0),
            "SVM": SVC(C=0.5, kernel="rbf", gamma="scale", class_weight=cw,
                       probability=True, random_state=RANDOM_STATE),
        }
        self.models = {n: self._make_pipeline(c) for n, c in base.items()}
        print(f"    {len(self.models)} pipelines: {list(self.models)}")
        return self.models

    # =============================================== 7. TRAIN + TEST + TABLES
    def train_and_evaluate(self):
        print("\n[7] Training and evaluating on the held-out test set ...")
        for name, pipe in self.models.items():
            print(f"    fitting {name} ...")
            pipe.fit(self.X_train, self.y_train)
            self.fitted[name] = pipe
            self.results[name] = dict(y_pred=pipe.predict(self.X_test),
                                      y_pred_proba=pipe.predict_proba(self.X_test))
        return self.results

    def calculate_metrics(self):
        print("\n[8] Performance of the cellitac models ...")
        rows = []
        for name, r in self.results.items():
            yp, pp = r["y_pred"], r["y_pred_proba"]
            try:
                auc = roc_auc_score(self.y_test, pp, multi_class="ovr",
                                    average="weighted")
            except Exception:
                auc = np.nan
            rows.append(dict(
                Model=name,
                Balanced_Acc=balanced_accuracy_score(self.y_test, yp),
                Macro_F1=f1_score(self.y_test, yp, average="macro"),
                Accuracy=accuracy_score(self.y_test, yp),
                Weighted_F1=f1_score(self.y_test, yp, average="weighted"),
                Precision=precision_score(self.y_test, yp, average="weighted"),
                Recall=recall_score(self.y_test, yp, average="weighted"),
                AUC=auc))
        self.metrics_df = pd.DataFrame(rows)
        print(self.metrics_df.round(4).to_string(index=False))
        self._save_table(self.metrics_df, "table02_model_performance",
                         "This table shows the performance of the cellitac models")
        return self.metrics_df

    def per_class_metrics(self):
        print("\n[9] Per-class metrics ...")
        frames = []
        for name, r in self.results.items():
            rep = classification_report(self.y_test, r["y_pred"],
                                        target_names=self.class_names,
                                        output_dict=True, zero_division=0)
            df = (pd.DataFrame(rep).T
                    .drop(index=["accuracy", "macro avg", "weighted avg"],
                          errors="ignore")
                    .reset_index().rename(columns={"index": "Cell_Type"}))
            df.insert(0, "Model", name)
            frames.append(df)
        self.per_class_df = pd.concat(frames, ignore_index=True)
        self.per_class_df["support"] = self.per_class_df["support"].astype(int)
        self._save_table(self.per_class_df, "table04_per_class_metrics",
                         "Per-class performance on the held-out test set")
        # compact recall matrix, one column per model
        recall = self.per_class_df.pivot(index="Cell_Type", columns="Model",
                                         values="recall").reset_index()
        self._save_table(recall, "table05_per_class_recall",
                         "Per-class recall by model (held-out test set)")
        return self.per_class_df

    # ====================================== 10. CV, OVERFITTING, ACC TABLE
    def detect_overfitting(self):
        print("\n[10] Cross-validation and overfitting score ...")
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
        rows = []
        for name, pipe in self.models.items():
            val_pred = cross_val_predict(pipe, self.X_train, self.y_train,
                                         cv=cv, n_jobs=-1)
            cv_bal = balanced_accuracy_score(self.y_train, val_pred)
            fitted = self.fitted[name]
            train_bal = balanced_accuracy_score(self.y_train,
                                                fitted.predict(self.X_train))
            gap = train_bal - cv_bal
            test_acc = accuracy_score(self.y_test, self.results[name]["y_pred"])
            verdict = ("Severe overfitting" if gap > 0.30 else
                       "High overfitting"   if gap > 0.20 else
                       "Moderate"           if gap > 0.10 else
                       "Mild"               if gap > 0.05 else "Good fit")
            rows.append(dict(Model=name,
                             Test_Accuracy_pct=round(test_acc * 100, 2),
                             Train_BalAcc=train_bal, CV_BalAcc=cv_bal,
                             Overfitting_Score=gap, Verdict=verdict))
        self.cv_results = pd.DataFrame(rows)
        print(self.cv_results.round(4).to_string(index=False))
        # exactly the table requested: Model | Test Accuracy (%) | Overfitting Score
        self._save_table(
            self.cv_results[["Model", "Test_Accuracy_pct",
                             "Overfitting_Score", "Verdict"]],
            "table03_accuracy_and_overfitting",
            "Test accuracy and overfitting score per model\n"
            "(overfitting score = train balanced accuracy - 5-fold CV balanced accuracy)")
        self._save_table(self.cv_results, "table03b_cv_full",
                         "Cross-validation detail (leakage-safe, pipeline re-fit per fold)")
        return self.cv_results

    # ============================================ 11. PER-MODEL TF IMPORTANCE
    def _permutation_importance_selected(self, name, feats):
        """Manual permutation importance restricted to the model's selected
        motifs. Used for estimators with no native importance (e.g. RBF-SVM)."""
        pipe = self.fitted[name]
        rng = np.random.default_rng(RANDOM_STATE)
        n = min(PERM_SUBSAMPLE, self.X_test.shape[0])
        idx = rng.choice(self.X_test.shape[0], size=n, replace=False)
        Xs, ys = self.X_test.iloc[idx], self.y_test[idx]
        base = balanced_accuracy_score(ys, pipe.predict(Xs))
        drops = {}
        for k, f in enumerate(feats, 1):
            vals = []
            for _ in range(PERM_N_REPEATS):
                Xp = Xs.copy()
                Xp[f] = rng.permutation(Xp[f].values)
                vals.append(base - balanced_accuracy_score(ys, pipe.predict(Xp)))
            drops[f] = float(np.mean(vals))
            if k % 25 == 0:
                print(f"        permutation {k}/{len(feats)} motifs")
        return pd.Series(drops).sort_values(ascending=False)

    def _model_importance(self, name):
        """Return (importance Series, method label) for one fitted model."""
        pipe = self.fitted[name]
        feats = self.X_train.columns[pipe.named_steps["select"].get_support()]
        clf = pipe.named_steps["clf"]
        if hasattr(clf, "feature_importances_"):
            return (pd.Series(clf.feature_importances_, index=feats)
                      .sort_values(ascending=False),
                    "native gain / impurity importance")
        if hasattr(clf, "coef_"):
            return (pd.Series(np.abs(clf.coef_).mean(axis=0), index=feats)
                      .sort_values(ascending=False),
                    "mean |coefficient| across classes")
        print(f"    {name}: no native importance -> permutation importance "
              f"({PERM_N_REPEATS} repeats, {min(PERM_SUBSAMPLE, len(self.y_test))} test cells)")
        return (self._permutation_importance_selected(name, list(feats)),
                "permutation importance (drop in balanced accuracy)")

    def feature_importance_per_model(self):
        """# NEW — top-20 bar chart with scores, for every model."""
        print("\n[11] TF importance per model ...")
        self.importance, self.importance_method = {}, {}
        for name in self.fitted:
            imp, method = self._model_importance(name)
            self.importance[name] = imp
            self.importance_method[name] = method
            top = imp.head(20)

            out = top.rename("Importance").reset_index()
            out.columns = ["TF_motif", "Importance"]
            out.insert(0, "Rank", np.arange(1, len(out) + 1))
            self._save_table(out, f"table06_top20_TF_{name.replace(' ', '_')}",
                             f"Top 20 TF motifs — {name}\n({method})",
                             float_fmt="{:.5f}")

            asc = top.sort_values(ascending=True)   # ascending -> largest on top
            fig, ax = plt.subplots(figsize=(10, 8))
            colors = sns.color_palette("viridis", len(asc))
            bars = ax.barh(range(len(asc)), asc.values, color=colors)
            ax.set_yticks(range(len(asc)))
            ax.set_yticklabels(asc.index, fontsize=10)
            span = float(asc.values.max()) if asc.values.max() > 0 else 1.0
            for b, v in zip(bars, asc.values):
                ax.text(b.get_width() + span * 0.012, b.get_y() + b.get_height() / 2,
                        f"{v:.4f}", va="center", fontsize=9)
            ax.set_xlim(0, span * 1.18)
            ax.set_xlabel(f"Importance — {method}")
            ax.set_title(f"Top 20 TF motifs — {name}\n"
                         f"features selected in-fold from "
                         f"{self.X_train.shape[1]} chromVAR motifs",
                         fontsize=13, weight="bold")
            ax.grid(axis="x", alpha=0.3)
            fig.tight_layout()
            self._save(fig, f"fig06_importance_{name.replace(' ', '_')}.png")

        combined = pd.concat(
            [s.head(20).rename("Importance").reset_index().assign(Model=n)
             for n, s in self.importance.items()], ignore_index=True)
        combined.columns = ["TF_motif", "Importance", "Model"]
        combined.to_csv(self.out_dir / "table06_top20_TF_all_models.csv", index=False)
        return self.importance

    # ================ 12. STATISTICAL TF <-> CELL-TYPE ASSOCIATION (# NEW)
    def tf_celltype_associations(self):
        """A real criterion, replacing the old 'mean + 1 SD' rule.

        For every TF motif x every cell type, on TRAIN cells only:
            * one-vs-rest Mann-Whitney U test  -> p value (no normality assumed)
            * rank-biserial AUC = U / (n1 * n2) -> scale-free specificity,
              0.5 = no separation, 1.0 = perfectly higher in that class
            * Cohen's d with pooled SD          -> standardised effect size
        p values are corrected across ALL TF x class tests with Benjamini-Hochberg.
        A TF is called characteristic of a class when
            FDR < FDR_ALPHA  AND  d >= MIN_COHENS_D  AND  AUC >= MIN_RB_AUC.
        """
        print("\n[12] Statistical TF <-> cell-type association (train cells only) ...")
        y_tr = pd.Series(self.label_encoder.inverse_transform(self.y_train),
                         index=self.X_train.index)
        recs = []
        for ct in self.class_names:
            m = (y_tr == ct).values
            n1, n2 = int(m.sum()), int((~m).sum())
            if n1 < 3 or n2 < 3:
                continue
            A, B = self.X_train.values[m], self.X_train.values[~m]
            mu1, mu2 = A.mean(0), B.mean(0)
            s1, s2 = A.std(0, ddof=1), B.std(0, ddof=1)
            pooled = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) /
                             (n1 + n2 - 2))
            pooled[pooled == 0] = np.nan
            d = (mu1 - mu2) / pooled
            U, p = stats.mannwhitneyu(A, B, alternative="two-sided", axis=0)
            auc = U / (n1 * n2)
            for j, tf in enumerate(self.X_train.columns):
                recs.append(dict(Cell_Type=ct, TF_motif=tf,
                                 Mean_in_class=float(mu1[j]),
                                 Mean_other=float(mu2[j]),
                                 Cohens_d=float(d[j]), RankBiserial_AUC=float(auc[j]),
                                 p_value=float(p[j]), n_in_class=n1))
        assoc = pd.DataFrame(recs)
        assoc["FDR"] = self._bh_fdr(assoc["p_value"].values)
        assoc["Significant"] = ((assoc["FDR"] < FDR_ALPHA) &
                                (assoc["Cohens_d"] >= MIN_COHENS_D) &
                                (assoc["RankBiserial_AUC"] >= MIN_RB_AUC))
        # association strength used for network edge length
        assoc["Assoc_Strength"] = np.where(assoc["Significant"], assoc["Cohens_d"], 0.0)
        assoc = assoc.sort_values(["Cell_Type", "Cohens_d"], ascending=[True, False])
        self.assoc = assoc

        sig = assoc[assoc["Significant"]]
        print(f"    tested {len(assoc):,} TF x class pairs | "
              f"{len(sig):,} significant at FDR<{FDR_ALPHA}, "
              f"d>={MIN_COHENS_D}, AUC>={MIN_RB_AUC}")
        print("    significant TFs per class:")
        for ct in self.class_names:
            k = int((sig["Cell_Type"] == ct).sum())
            top = sig[sig["Cell_Type"] == ct].head(5)["TF_motif"].tolist()
            print(f"      {ct:<20} {k:>4}   top: {', '.join(top) if top else '-'}")

        assoc.to_csv(self.out_dir / "table07_TF_celltype_associations_full.csv",
                     index=False)
        top_sig = (sig.groupby("Cell_Type", group_keys=False)
                      .head(10)
                      .loc[:, ["Cell_Type", "TF_motif", "Cohens_d",
                               "RankBiserial_AUC", "FDR", "n_in_class"]]
                      .reset_index(drop=True))
        self._save_table(top_sig, "table07_TF_celltype_associations_top",
                         "TF motifs characteristic of each cell type\n"
                         f"(Mann-Whitney U one-vs-rest, BH-FDR < {FDR_ALPHA}, "
                         f"Cohen's d >= {MIN_COHENS_D}, rank-biserial AUC >= {MIN_RB_AUC})")
        return assoc

    # ============================================ 13. TF -> CELL-TYPE NETWORK
    def build_networks(self):
        """# NEW — one network per model. Edge length encodes association
        strength: the stronger the association, the shorter the edge."""
        print("\n[13] TF -> cell-type networks (one per model) ...")
        self.networks = {}
        for name in self.fitted:
            cand = list(self.importance[name].head(NET_TOP_TFS).index)
            sub = self.assoc[(self.assoc["Significant"]) &
                             (self.assoc["TF_motif"].isin(cand))]
            if sub.empty:
                print(f"    [skip] {name}: no significant TF-cell association "
                      "among its top motifs")
                continue

            G = nx.Graph()
            used_tfs = sorted(sub["TF_motif"].unique())
            for tf in used_tfs:
                G.add_node(tf, node_type="TF")
            for ct in sorted(sub["Cell_Type"].unique()):
                G.add_node(ct, node_type="cell_type")
            smax = float(sub["Assoc_Strength"].max())
            for _, r in sub.iterrows():
                G.add_edge(r["TF_motif"], r["Cell_Type"], kind="tf-ct",
                           strength=float(r["Assoc_Strength"]),
                           weight=float(r["Assoc_Strength"]) / smax,
                           auc=float(r["RankBiserial_AUC"]),
                           fdr=float(r["FDR"]))
            if SHOW_TF_TF_EDGES and len(used_tfs) > 1:
                cm = self.X_train[used_tfs].corr()
                for i, a in enumerate(used_tfs):
                    for b in used_tfs[i + 1:]:
                        if abs(cm.loc[a, b]) > 0.5:
                            G.add_edge(a, b, kind="tf-tf",
                                       strength=abs(cm.loc[a, b]),
                                       weight=abs(cm.loc[a, b]) * 0.3)

            # spring layout: higher weight -> stronger pull -> shorter edge
            pos = nx.spring_layout(G, weight="weight", k=1.3, iterations=800,
                                   seed=RANDOM_STATE, scale=1.0)
            pos = self._resolve_overlaps(pos, min_dist=0.36)
            fig, ax = plt.subplots(figsize=(17, 12))
            tf_nodes = [n for n, d in G.nodes(data=True) if d["node_type"] == "TF"]
            ct_nodes = [n for n, d in G.nodes(data=True) if d["node_type"] == "cell_type"]

            tf_ct = [(u, v, d) for u, v, d in G.edges(data=True) if d["kind"] == "tf-ct"]
            if SHOW_TF_TF_EDGES:
                tf_tf = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == "tf-tf"]
                nx.draw_networkx_edges(G, pos, edgelist=tf_tf, style="dashed",
                                       edge_color="#9aa0a6", alpha=0.35, width=0.8, ax=ax)
            strengths = np.array([d["strength"] for _, _, d in tf_ct])
            widths = 1.0 + 4.0 * (strengths - strengths.min()) / \
                     (np.ptp(strengths) + 1e-9)
            edges = nx.draw_networkx_edges(
                G, pos, edgelist=[(u, v) for u, v, _ in tf_ct], width=widths,
                edge_color=strengths, edge_cmap=plt.cm.YlOrRd,
                edge_vmin=float(strengths.min()), edge_vmax=float(strengths.max()),
                alpha=0.9, ax=ax)

            deg = dict(G.degree())
            nx.draw_networkx_nodes(G, pos, nodelist=tf_nodes, node_color="#7FB3D5",
                                   node_size=[420 + 130 * deg[n] for n in tf_nodes],
                                   edgecolors="#2E4A62", linewidths=1.2, ax=ax)
            nx.draw_networkx_nodes(G, pos, nodelist=ct_nodes, node_color="#E8746A",
                                   node_size=[1100 + 160 * deg[n] for n in ct_nodes],
                                   edgecolors="#8C2F26", linewidths=1.6, ax=ax)
            span = max(np.ptp([p[1] for p in pos.values()]), 1e-6)
            for n in tf_nodes:                      # TF labels sit below the node
                ax.text(pos[n][0], pos[n][1] - 0.052 * span, n, fontsize=8,
                        ha="center", va="top", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                  ec="#7FB3D5", lw=0.8, alpha=0.9))
            for n in ct_nodes:                      # cell-type labels above
                ax.text(pos[n][0], pos[n][1] + 0.075 * span, n, fontsize=10.5,
                        weight="bold", ha="center", va="bottom", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.2", fc="#FDEDEB",
                                  ec="#8C2F26", lw=1.2, alpha=0.95))
            ax.margins(0.16)
            ax.set_aspect("equal", adjustable="datalim")

            cbar = fig.colorbar(edges, ax=ax, fraction=0.035, pad=0.02)
            cbar.set_label("Association strength (Cohen's d)", fontsize=10)

            from matplotlib.lines import Line2D
            from matplotlib.patches import Patch
            key = [
                Line2D([], [], marker="o", color="none", markerfacecolor="#7FB3D5",
                       markeredgecolor="#2E4A62", markersize=12,
                       label="TF motif (size = number of cell types linked)"),
                Line2D([], [], marker="o", color="none", markerfacecolor="#E8746A",
                       markeredgecolor="#8C2F26", markersize=15,
                       label="Cell type (size = number of TFs linked)"),
                Line2D([], [], color="#B22222", lw=4,
                       label="Edge: significant TF-cell association"),
                Line2D([], [], color="#FDBE85", lw=1.2,
                       label="Thin / pale edge = weaker association"),
                Patch(facecolor="white", edgecolor="white",
                      label="Edge LENGTH is inverse to strength:\n"
                            "shorter edge = the TF characterises that\n"
                            "cell type more strongly"),
                Patch(facecolor="white", edgecolor="white",
                      label=f"Criterion: BH-FDR < {FDR_ALPHA}, Cohen's d >= "
                            f"{MIN_COHENS_D},\nrank-biserial AUC >= {MIN_RB_AUC}"),
            ]
            ax.legend(handles=key, loc="upper left", bbox_to_anchor=(1.06, 1.0),
                      fontsize=9, frameon=True, title="Key", title_fontsize=11)
            ax.set_title(
                f"TF motif -> cell-type regulatory network — {name}\n"
                f"{len(tf_nodes)} TF motifs, {len(ct_nodes)} cell types, "
                f"{len(tf_ct)} significant associations "
                f"(candidates = top {NET_TOP_TFS} motifs of this model)",
                fontsize=13, weight="bold")
            ax.axis("off")
            fig.tight_layout()
            self._save(fig, f"fig07_TF_network_{name.replace(' ', '_')}.png")

            edge_tbl = pd.DataFrame([
                dict(Model=name, TF_motif=u if G.nodes[u]["node_type"] == "TF" else v,
                     Cell_Type=v if G.nodes[v]["node_type"] == "cell_type" else u,
                     Cohens_d=d["strength"], RankBiserial_AUC=d["auc"], FDR=d["fdr"])
                for u, v, d in tf_ct])
            edge_tbl.to_csv(
                self.out_dir / f"table08_network_edges_{name.replace(' ', '_')}.csv",
                index=False)
            self.networks[name] = G
            print(f"    {name}: {G.number_of_nodes()} nodes, "
                  f"{G.number_of_edges()} edges")
        return self.networks

    # ==================================================== 14. LEARNING CURVES
    def plot_learning_curves(self):
        """# NEW — real learning curves from the leakage-safe pipeline."""
        print("\n[14] Learning curves (balanced accuracy, pipeline re-fit per fold) ...")
        cv = StratifiedKFold(n_splits=LC_CV_FOLDS, shuffle=True,
                             random_state=RANDOM_STATE)
        n_models = len(self.models)
        fig, axes = plt.subplots(1, n_models, figsize=(5.2 * n_models, 4.8),
                                 sharey=True)
        axes = np.atleast_1d(axes)
        rows = []
        for ax, (name, pipe) in zip(axes, self.models.items()):
            print(f"    {name} ...")
            sizes, train_sc, val_sc = learning_curve(
                pipe, self.X_train, self.y_train, cv=cv,
                train_sizes=LC_TRAIN_SIZES, scoring="balanced_accuracy",
                n_jobs=-1, shuffle=True, random_state=RANDOM_STATE)
            tm, ts = train_sc.mean(1), train_sc.std(1)
            vm, vs = val_sc.mean(1), val_sc.std(1)
            ax.plot(sizes, tm, "o-", color="#C0392B", label="Training")
            ax.fill_between(sizes, tm - ts, tm + ts, alpha=0.15, color="#C0392B")
            ax.plot(sizes, vm, "s-", color="#2874A6", label="Cross-validation")
            ax.fill_between(sizes, vm - vs, vm + vs, alpha=0.15, color="#2874A6")
            ax.set_title(f"{name}\nfinal gap = {tm[-1] - vm[-1]:.3f}",
                         fontsize=11, weight="bold")
            ax.set_xlabel("Training cells")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=9, loc="lower right")
            for s, a, b in zip(sizes, tm, vm):
                rows.append(dict(Model=name, Train_Size=int(s),
                                 Train_BalAcc=a, CV_BalAcc=b, Gap=a - b))
        axes[0].set_ylabel("Balanced accuracy")
        fig.suptitle("Learning curves — cellitac (chromVAR TF motif activity)",
                     fontsize=14, weight="bold")
        fig.tight_layout()
        self._save(fig, "fig08_learning_curves.png")
        self.lc_df = pd.DataFrame(rows)
        self.lc_df.to_csv(self.out_dir / "table09_learning_curve.csv", index=False)
        return self.lc_df

    # ============================== 15. UMAP + t-SNE OF THE LEARNED SPACE
    def plot_learned_embedding(self, model=None):
        """# NEW — after training: embed the TEST cells in the model's own
        selected-feature space and colour by true vs predicted label."""
        if model is None:
            model = self.metrics_df.loc[self.metrics_df["Macro_F1"].idxmax(), "Model"]
        print(f"\n[15] UMAP + t-SNE of the learned space ({model}) ...")
        pipe = self.fitted[model]
        Xs = pipe.named_steps["scale"].transform(self.X_test)
        Xs = pipe.named_steps["select"].transform(Xs)     # model's own subspace
        y_true = self.label_encoder.inverse_transform(self.y_test)
        y_pred = self.label_encoder.inverse_transform(self.results[model]["y_pred"])

        embs = {}
        try:
            import umap
            embs["UMAP"] = umap.UMAP(n_neighbors=30, min_dist=0.3,
                                     random_state=RANDOM_STATE).fit_transform(Xs)
        except Exception as e:
            print(f"    [warn] umap unavailable ({e}); UMAP panel skipped")
        embs["t-SNE"] = TSNE(n_components=2, perplexity=30, init="pca",
                             random_state=RANDOM_STATE).fit_transform(Xs)

        colors = sns.color_palette(PALETTE, len(self.class_names))
        cmap = {ct: colors[i] for i, ct in enumerate(self.class_names)}
        fig, axes = plt.subplots(len(embs), 2, figsize=(15, 6.6 * len(embs)),
                                 squeeze=False)
        correct = (y_true == y_pred)
        for r, (mname, emb) in enumerate(embs.items()):
            for c, (lab, ttl) in enumerate([(y_true, "True label"),
                                            (y_pred, "Predicted label")]):
                ax = axes[r][c]
                for ct in self.class_names:
                    m = lab == ct
                    if m.sum():
                        ax.scatter(emb[m, 0], emb[m, 1], s=7, alpha=0.7,
                                   color=cmap[ct], linewidths=0, label=ct)
                ax.set_title(f"{mname} — {ttl}", fontsize=12, weight="bold")
                ax.set_xlabel(f"{mname} 1"); ax.set_ylabel(f"{mname} 2")
                if r == 0 and c == 1:
                    ax.legend(markerscale=3, fontsize=8, loc="center left",
                              bbox_to_anchor=(1.01, 0.5), frameon=False)
        fig.suptitle(
            f"Separation achieved by {model} in its own selected-feature space\n"
            f"{Xs.shape[0]:,} held-out test cells x {Xs.shape[1]} selected TF motifs "
            f"— {correct.mean() * 100:.1f}% of cells correctly assigned",
            fontsize=14, weight="bold")
        fig.tight_layout()
        self._save(fig, "fig09_umap_tsne_after_training.png")
        return embs

    # ======================================== 16. STANDARD COMPARISON PLOTS
    def comparison_plots(self):
        print("\n[16] Model comparison and confusion matrices ...")
        m = self.metrics_df.melt(id_vars="Model", var_name="Metric",
                                 value_name="Score")
        fig, ax = plt.subplots(figsize=(11, 6))
        sns.barplot(data=m[m.Metric.isin(["Balanced_Acc", "Macro_F1", "AUC"])],
                    x="Model", y="Score", hue="Metric", ax=ax)
        for cont in ax.containers:
            ax.bar_label(cont, fmt="%.3f", fontsize=8, padding=2)
        ax.set_ylim(0, 1.05)
        ax.set_title("Model comparison (imbalance-aware metrics)",
                     fontsize=13, weight="bold")
        plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
        fig.tight_layout()
        self._save(fig, "fig04_model_comparison.png")

        n = len(self.results)
        fig, axes = plt.subplots(1, n, figsize=(5.6 * n, 5.2))
        axes = np.atleast_1d(axes)
        for ax, (name, r) in zip(axes, self.results.items()):
            cm = confusion_matrix(self.y_test, r["y_pred"])
            sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False,
                        xticklabels=self.class_names, yticklabels=self.class_names)
            ax.set_title(name, fontsize=12, weight="bold")
            ax.set_xlabel("Predicted"); ax.set_ylabel("True")
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
            plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
        fig.suptitle("Confusion matrices on the held-out test set",
                     fontsize=14, weight="bold")
        fig.tight_layout()
        self._save(fig, "fig05_confusion_matrices.png")
        return True

    # ========================================= 17. FEATURE-SET BASELINES
    def compare_feature_baselines(self):
        print("\n[17] Baseline comparison (chromVAR vs RNA-only vs peak LSI) ...")
        files = {"ATAC_chromVAR_TF": "cellitac_TF_activity.csv",
                 "RNA_only":         "rna_baseline_features.csv",
                 "ATAC_peak_LSI":    "atac_peak_baseline_lsi.csv"}
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
        y = pd.Series(self.y_encoded, index=self.X.index)
        rows = []
        for label, fname in files.items():
            f = self.data_dir / fname
            if not f.exists():
                print(f"    [skip] {fname} not found")
                continue
            Xb = pd.read_csv(f, index_col=0)
            cells = Xb.index.intersection(self.X.index)
            Xb, yb = Xb.loc[cells], y.loc[cells].values
            steps = [("scale", StandardScaler()),
                     ("select", SelectKBest(f_classif,
                                            k=min(MAX_FEATURES, Xb.shape[1])))]
            if USE_SMOTE:
                steps.append(("smote", SMOTE(k_neighbors=3, random_state=RANDOM_STATE)))
            steps.append(("clf", LogisticRegression(
                penalty="l2", C=1.0,
                class_weight=(None if USE_SMOTE else "balanced"),
                max_iter=3000, random_state=RANDOM_STATE)))
            pred = cross_val_predict(ImbPipeline(steps), Xb, yb, cv=cv, n_jobs=-1)
            rows.append(dict(Feature_Set=label, n_features=int(Xb.shape[1]),
                             CV_BalAcc=balanced_accuracy_score(yb, pred),
                             CV_Macro_F1=f1_score(yb, pred, average="macro")))
        self.baseline_df = pd.DataFrame(rows)
        print(self.baseline_df.round(4).to_string(index=False))
        self._save_table(self.baseline_df, "table10_feature_set_baselines",
                         "Same pipeline, same cells, same folds — three feature sets")
        return self.baseline_df

    # ==================================================== 18. FINAL REPORT
    def final_report(self):
        print("\n[18] Final report ...")
        best = self.metrics_df.loc[self.metrics_df["Macro_F1"].idxmax()]
        sig = self.assoc[self.assoc["Significant"]]
        report = dict(
            dataset=dict(cells_after_QC=int(self.counts_before.sum()),
                         cells_modelled=int(len(self.y)),
                         tf_motifs_input=int(self.X.shape[1]),
                         classes=self.class_names,
                         class_counts={k: int(v) for k, v in
                                       self.counts_after.items()}),
            balancing="SMOTE(in-fold)" if USE_SMOTE else "class_weight=balanced",
            feature_selection=dict(max_features=MAX_FEATURES,
                                   fitted="inside each training fold"),
            performance=self.metrics_df.round(4).to_dict("records"),
            cross_validation=self.cv_results.round(4).to_dict("records"),
            tf_association_criterion=dict(test="Mann-Whitney U, one-vs-rest",
                                          correction="Benjamini-Hochberg",
                                          fdr_alpha=FDR_ALPHA,
                                          min_cohens_d=MIN_COHENS_D,
                                          min_rank_biserial_auc=MIN_RB_AUC,
                                          n_tests=int(len(self.assoc)),
                                          n_significant=int(len(sig))),
            importance_method={k: v for k, v in self.importance_method.items()},
            best_model=dict(name=best["Model"],
                            macro_f1=float(best["Macro_F1"]),
                            balanced_acc=float(best["Balanced_Acc"])),
        )
        with open(self.out_dir / "cellitac_ml_report.json", "w") as f:
            json.dump(report, f, indent=2, default=str)
        print(f"    best model by macro-F1: {best['Model']} "
              f"(macro-F1 {best['Macro_F1']:.3f}, bal-acc {best['Balanced_Acc']:.3f})")
        return report

    # ============================================================ RUN ALL
    def run(self):
        steps = [
            self.load_data, self.analyze_class_imbalance, self.clean_features,
            self.plot_input_embedding, self.prepare_train_test_split,
            self.define_models, self.train_and_evaluate, self.calculate_metrics,
            self.per_class_metrics, self.detect_overfitting,
            self.feature_importance_per_model, self.tf_celltype_associations,
            self.build_networks, self.plot_learning_curves,
            self.plot_learned_embedding, self.comparison_plots,
            self.compare_feature_baselines, self.final_report,
        ]
        for i, step in enumerate(steps, 1):
            try:
                step()
            except Exception as e:
                print(f"\n[FAIL] step {i}/{len(steps)} ({step.__name__}): {e}")
                import traceback
                traceback.print_exc()
                return False
        print("\n" + "=" * 78)
        print("cellitac v3 completed. Figures and tables written to:",
              self.out_dir.resolve())
        print("=" * 78)
        return True
def main():
    ap = argparse.ArgumentParser(description="cellitac TF-centric ML pipeline v3")
    
    # هنا عدلنا الـ default لكل واحد فيهم للمسارات الجديدة
    ap.add_argument(
        "--data-dir", 
        default="cellitac_output", 
        help="folder with the R outputs"
    )
    ap.add_argument(
        "--out-dir", 
        default="output/output_Second_Framework_after_dropping", 
        help="output folder"
    )
    
    args = ap.parse_args()
    pipe = CellitacPipeline(data_dir=args.data_dir, out_dir=args.out_dir)
    
    if pipe.run():
        print("\nHeadline results:")
        print(pipe.metrics_df.round(4).to_string(index=False))
        print()
        print(pipe.cv_results[["Model", "Test_Accuracy_pct",
                               "Overfitting_Score", "Verdict"]].to_string(index=False))

if __name__ == "__main__":
    main()