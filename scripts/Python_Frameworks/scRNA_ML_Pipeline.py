#!/usr/bin/env python3
"""
Complete Machine Learning Pipeline for scRNA-seq Cell Type Classification
RNA-ONLY VERSION - For comparison with integrated RNA+ATAC results
Handles class imbalance, prevents overfitting, and provides comprehensive evaluation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import openpyxl
import json
from pathlib import Path
import numpy as np
from math import pi
import networkx as nx

import warnings
warnings.filterwarnings('ignore')

# ML Libraries
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier

# Preprocessing Libraries
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif, RFE
from sklearn.model_selection import learning_curve
from sklearn.decomposition import PCA

# Imbalance handling Libraries
from imblearn.over_sampling import SMOTE, BorderlineSMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.combine import SMOTETomek

# Evaluation Libraries
from sklearn.metrics import (classification_report, confusion_matrix, 
                           accuracy_score, precision_score, recall_score, 
                           f1_score, roc_auc_score, roc_curve, auc)

# Visualization Libraries
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# Define mint-green color for RNA-only results
MINT_GREEN = '#90EE90'
MINT_GREEN_DARK = '#3CB371'

class scRNAMLPipeline:
    """Complete ML Pipeline for scRNA-seq ONLY cell type classification"""
    
    def __init__(self, data_dir):
        """Initialize pipeline with data directory"""
        self.data_dir = Path(data_dir)
        self.models = {}
        self.results = {}
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        print("=== scRNA-seq ML Pipeline Initialized (RNA-ONLY) ===")

    def load_data(self):
        """Load and prepare RNA-ONLY data"""
        print("\n1. Loading RNA-seq data only...")
        
        try:
            # Load ONLY RNA features
            rna_file = self.data_dir / "rna_features_2000.csv"
            if not rna_file.exists():
                raise FileNotFoundError(f"RNA features file not found: {rna_file}")
            
            self.X = pd.read_csv(rna_file, index_col=0)
            print(f"   RNA features loaded: {self.X.shape}")
            print(f"   Data type: scRNA-seq ONLY")
            
            # Load labels
            if (self.data_dir / "cell_labels.csv").exists():
                labels_df = pd.read_csv(self.data_dir / "cell_labels.csv", index_col=0)
            else:
                labels_df = pd.read_csv(self.data_dir / "cell_type_distribution.csv", index_col=0)
            
            # Match labels to features
            common_cells = self.X.index.intersection(labels_df.index)
            self.X = self.X.loc[common_cells]
            self.y = labels_df.loc[common_cells, 'cell_type']
            
            print(f"   Final dataset: {self.X.shape[0]} cells, {self.X.shape[1]} features (RNA genes)")
            print(f"   Cell types: {self.y.value_counts().to_dict()}")

            # Remove rare classes (<10 samples)
            counts = self.y.value_counts()
            keep_types = counts[counts >= 10].index
            mask = self.y.isin(keep_types)
            removed = (~mask).sum()
            
            if removed > 0:
                rare_types = counts[counts < 10].index.tolist()
                self.X = self.X.loc[mask]
                self.y = self.y.loc[mask]
                print(f"   Removed {removed} cells from rare types: {rare_types}")
            else:
                print("   No rare cell types (<10) to remove.")
                
            if self.X.empty:
                raise ValueError("All cells were filtered out by the rare-type threshold.")

            # Encode labels for the model
            self.y_encoded = self.label_encoder.fit_transform(self.y)
            self.class_names = self.label_encoder.classes_
            return True
            
        except Exception as e:
            print(f"   Error loading data: {e}")
            return False

    def analyze_class_imbalance(self):
        """Analyze and visualize class distribution"""
        print("\n2. Analyzing class distribution...")
        
        # Calculate class distribution
        class_counts = pd.Series(self.y).value_counts()
        class_props = pd.Series(self.y).value_counts(normalize=True)
        
        # Create imbalance analysis
        imbalance_df = pd.DataFrame({
            'Cell_Type': class_counts.index,
            'Count': class_counts.values,
            'Proportion': class_props.values,
            'Imbalance_Ratio': class_counts.max() / class_counts.values
        })
        
        print("   Class Distribution:")
        print(imbalance_df.round(3))
        
        # Check if balancing is needed
        max_ratio = imbalance_df['Imbalance_Ratio'].max()
        self.needs_balancing = max_ratio > 2.0
        
        print(f"   Maximum imbalance ratio: {max_ratio:.2f}")
        print(f"   Balancing needed: {self.needs_balancing}")
        
        # Visualizations
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Class counts
        sns.barplot(data=imbalance_df, x='Cell_Type', y='Count', ax=axes[0,0], color=MINT_GREEN)
        axes[0,0].set_title('Cell Type Distribution (Counts) - RNA-only')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Class proportions
        sns.barplot(data=imbalance_df, x='Cell_Type', y='Proportion', ax=axes[0,1], color=MINT_GREEN)
        axes[0,1].set_title('Cell Type Distribution (Proportions) - RNA-only')
        axes[0,1].tick_params(axis='x', rotation=45)
        
        # Imbalance ratios
        sns.barplot(data=imbalance_df, x='Cell_Type', y='Imbalance_Ratio', ax=axes[1,0], color=MINT_GREEN)
        axes[1,0].set_title('Class Imbalance Ratios - RNA-only')
        axes[1,0].axhline(y=2.0, color='red', linestyle='--', label='Threshold (2.0)')
        axes[1,0].tick_params(axis='x', rotation=45)
        axes[1,0].legend()
        
        # Pie chart
        axes[1,1].pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%', colors=plt.cm.Greens(np.linspace(0.4, 0.8, len(class_counts))))
        axes[1,1].set_title('Cell Type Distribution (Pie Chart) - RNA-only')
        
        plt.tight_layout()
        plt.savefig('class_distribution_analysis_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        self.imbalance_df = imbalance_df
        return imbalance_df
    
    def feature_engineering(self):
        """Apply feature engineering and selection"""
        print("\n3. Feature engineering...")
        
        # Remove features with zero variance
        zero_var_features = self.X.columns[self.X.var() == 0]
        if len(zero_var_features) > 0:
            self.X = self.X.drop(columns=zero_var_features)
            print(f"   Removed {len(zero_var_features)} zero-variance features")
        
        # Remove highly correlated features
        corr_matrix = self.X.corr().abs()
        upper_triangle = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        high_corr_features = [column for column in upper_triangle.columns if any(upper_triangle[column] > 0.95)]
        
        if len(high_corr_features) > 0:
            self.X = self.X.drop(columns=high_corr_features)
            print(f"   Removed {len(high_corr_features)} highly correlated features (>0.95)")
        
        # Feature selection using SelectKBest
        selector = SelectKBest(score_func=f_classif, k=min(1000, self.X.shape[1]//2))
        X_selected = selector.fit_transform(self.X, self.y_encoded)
        selected_features = self.X.columns[selector.get_support()]
        
        self.X_engineered = pd.DataFrame(X_selected, index=self.X.index, columns=selected_features)
        
        print(f"   Selected top {len(selected_features)} features using univariate selection")
        print(f"   Final feature matrix: {self.X_engineered.shape}")
        
        return self.X_engineered
    
    def prepare_train_test_split(self):
        """Create train/test splits with proper stratification"""
        print("\n4. Creating train/test splits...")
        
        # Stratified split maintaining class proportions
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            self.X_engineered, self.y_encoded,
            test_size=0.2, 
            random_state=42,
            stratify=self.y_encoded
        )
        
        print(f"   Training set: {self.X_train.shape[0]} samples")
        print(f"   Test set: {self.X_test.shape[0]} samples")
        
        # Show class distribution in splits
        train_dist = pd.Series(self.y_train).value_counts()
        test_dist = pd.Series(self.y_test).value_counts()
        
        split_df = pd.DataFrame({
            'Train_Count': train_dist,
            'Test_Count': test_dist,
            'Train_Prop': train_dist / len(self.y_train),
            'Test_Prop': test_dist / len(self.y_test)
        }).fillna(0)
        
        print("   Class distribution in splits:")
        print(split_df.round(3))
        
        return True
    
    def handle_class_imbalance(self):
        """Apply class balancing techniques if needed"""
        print("\n5. Handling class imbalance...")
        
        if not self.needs_balancing:
            print("   No significant imbalance detected. Skipping balancing.")
            self.X_train_balanced = self.X_train.copy()
            self.y_train_balanced = self.y_train.copy()
            return
        
        print("   Applying SMOTE for class balancing...")
        
        # Apply SMOTE
        smote = BorderlineSMOTE(kind="borderline-1", k_neighbors=3)
        self.X_train_balanced, self.y_train_balanced = smote.fit_resample(
            self.X_train, self.y_train
        )
        
        # Compare before and after
        before_counts = pd.Series(self.y_train).value_counts()
        after_counts = pd.Series(self.y_train_balanced).value_counts()
        
        comparison_df = pd.DataFrame({
            'Before_SMOTE': before_counts,
            'After_SMOTE': after_counts,
            'Improvement': after_counts / before_counts
        }).fillna(0)
        
        print("   Balancing results:")
        print(comparison_df)
        
        # Visualize balancing effect
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Before balancing
        before_counts.plot(kind='bar', ax=ax1, color='lightcoral')
        ax1.set_title('Class Distribution Before SMOTE - RNA-only')
        ax1.set_ylabel('Count')
        ax1.tick_params(axis='x', rotation=45)
        
        # After balancing
        after_counts.plot(kind='bar', ax=ax2, color=MINT_GREEN)
        ax2.set_title('Class Distribution After SMOTE - RNA-only')
        ax2.set_ylabel('Count')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig('class_balancing_comparison_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        return comparison_df
    
    def scale_features(self):
        """Scale features for algorithms that require it"""
        print("\n6. Scaling features...")
        
        # Fit scaler on training data only
        self.X_train_scaled = pd.DataFrame(
            self.scaler.fit_transform(self.X_train_balanced),
            index=self.X_train_balanced.index,
            columns=self.X_train_balanced.columns
        )
        
        self.X_test_scaled = pd.DataFrame(
            self.scaler.transform(self.X_test),
            index=self.X_test.index,
            columns=self.X_test.columns
        )
        
        print(f"   Features scaled using StandardScaler")
        print(f"   Training mean: {self.X_train_scaled.mean().mean():.6f}")
        print(f"   Training std: {self.X_train_scaled.std().mean():.6f}")
        
        return True

    def define_models(self):
        """Define models with optimized hyperparameters for preventing overfitting"""
        print("\n7. Defining ML models...")
        
        n_classes = len(self.class_names)
        n_features = self.X_train_balanced.shape[1]
        n_samples = self.X_train_balanced.shape[0]
        
        # Random Forest - Conservative parameters
        self.models['Random Forest'] = {
            'model': RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                min_samples_split=20,
                min_samples_leaf=10,
                max_features='sqrt',
                class_weight='balanced',
                random_state=42,
                n_jobs=-1
            ),
            'requires_scaling': False
        }
        
        # XGBoost - With regularization
        self.models['XGBoost'] = {
            'model': XGBClassifier(
                n_estimators=100,
                max_depth=6,
                learning_rate=0.1,
                subsample=0.8,
                colsample_bytree=0.8,
                reg_alpha=0.1,
                reg_lambda=1.0,
                min_child_weight=3,
                random_state=42,
                eval_metric='mlogloss' if n_classes > 2 else 'logloss',
                verbosity=0
            ),
            'requires_scaling': False
        }
        
        # SVM - With regularization
        self.models['SVM'] = {
            'model': SVC(
                C=1.0,
                kernel='rbf',
                gamma='scale',
                class_weight='balanced',
                probability=True,
                random_state=42
            ),
            'requires_scaling': True
        }
        
        print(f"   Defined {len(self.models)} models")
        for name, config in self.models.items():
            print(f"     - {name}: {type(config['model']).__name__}")
        
        return self.models
    
    def train_and_evaluate_models(self):
        """Train all models and evaluate performance"""
        print("\n8. Training and evaluating models...")
        
        for model_name, config in self.models.items():
            print(f"\n   Training {model_name}...")
            
            # Select appropriate data (scaled vs unscaled)
            if config['requires_scaling']:
                X_train = self.X_train_scaled
                X_test = self.X_test_scaled
            else:
                X_train = self.X_train_balanced
                X_test = self.X_test
            
            # Train model
            model = config['model']
            model.fit(X_train, self.y_train_balanced)
            
            # Predictions
            y_pred = model.predict(X_test)
            y_pred_proba = model.predict_proba(X_test)
            
            # Store results
            self.results[model_name] = {
                'model': model,
                'y_pred': y_pred,
                'y_pred_proba': y_pred_proba,
                'X_test': X_test
            }
            
            print(f"     {model_name} training completed")
        
        return self.results
    
    def calculate_metrics(self):
        """Calculate comprehensive evaluation metrics"""
        print("\n9. Calculating evaluation metrics...")
        
        metrics_data = []
        
        for model_name, results in self.results.items():
            y_pred = results['y_pred']
            y_pred_proba = results['y_pred_proba']
            
            # Basic metrics
            accuracy = accuracy_score(self.y_test, y_pred)
            precision = precision_score(self.y_test, y_pred, average='weighted')
            recall = recall_score(self.y_test, y_pred, average='weighted')
            f1 = f1_score(self.y_test, y_pred, average='weighted')
            
            # Multi-class AUC
            if len(self.class_names) > 2:
                auc_score = roc_auc_score(self.y_test, y_pred_proba, multi_class='ovr', average='weighted')
            else:
                auc_score = roc_auc_score(self.y_test, y_pred_proba[:, 1])
            
            metrics_data.append({
                'Model': model_name,
                'Accuracy': accuracy,
                'Precision': precision,
                'Recall': recall,
                'F1_Score': f1,
                'AUC': auc_score
            })
        
        self.metrics_df = pd.DataFrame(metrics_data)
        
        print("\n   Model Performance Summary (RNA-only):")
        print(self.metrics_df.round(4))
        
        return self.metrics_df
    
    def create_detailed_evaluation_table(self):
        """Create detailed per-class evaluation metrics"""
        print("\n10. Creating detailed evaluation tables...")
        
        detailed_results = {}
        
        for model_name, results in self.results.items():
            y_pred = results['y_pred']
            
            # Get per-class metrics
            report = classification_report(
                self.y_test, y_pred,
                target_names=self.class_names,
                output_dict=True
            )
            
            # Convert to DataFrame
            class_metrics = pd.DataFrame(report).transpose()
            class_metrics = class_metrics.drop(['accuracy'], errors='ignore')
            class_metrics = class_metrics.round(4)
            
            detailed_results[model_name] = class_metrics
            
            print(f"\n   {model_name} - Per-class Metrics:")
            print(class_metrics)
        
        self.detailed_results = detailed_results
        return detailed_results
    
    def create_visualizations(self):
        """Create comprehensive visualizations"""
        print("\n11. Creating visualizations...")
        
        # 1. Model comparison
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.ravel()
        
        metrics = ['Accuracy', 'Precision', 'Recall', 'F1_Score', 'AUC']
        
        for i, metric in enumerate(metrics):
            sns.barplot(data=self.metrics_df, x='Model', y=metric, ax=axes[i], color=MINT_GREEN)
            axes[i].set_title(f'{metric} Comparison - RNA-only')
            axes[i].set_ylim(0, 1)
            axes[i].tick_params(axis='x', rotation=45)
            
            # Add value labels
            for j, v in enumerate(self.metrics_df[metric]):
                axes[i].text(j, v + 0.01, f'{v:.3f}', ha='center', va='bottom')
        
        # Overall comparison
        metrics_melted = self.metrics_df.melt(id_vars=['Model'], var_name='Metric', value_name='Score')
        sns.barplot(data=metrics_melted, x='Model', y='Score', hue='Metric', ax=axes[5])
        axes[5].set_title('Overall Model Comparison - RNA-only')
        axes[5].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        axes[5].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig('model_performance_comparison_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # 2. Confusion matrices
        fig, axes = plt.subplots(1, len(self.models), figsize=(6*len(self.models), 5))
        if len(self.models) == 1:
            axes = [axes]
        
        for i, (model_name, results) in enumerate(self.results.items()):
            cm = confusion_matrix(self.y_test, results['y_pred'])
            sns.heatmap(cm, annot=True, fmt='d', cmap='Greens', 
                       xticklabels=self.class_names, 
                       yticklabels=self.class_names, 
                       ax=axes[i])
            axes[i].set_title(f'{model_name}\nConfusion Matrix - RNA-only')
            axes[i].set_ylabel('True Label')
            axes[i].set_xlabel('Predicted Label')
        
        plt.tight_layout()
        plt.savefig('confusion_matrices_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # 3. ROC Curves (if binary classification)
        if len(self.class_names) == 2:
            plt.figure(figsize=(10, 8))
            
            for model_name, results in self.results.items():
                y_pred_proba = results['y_pred_proba'][:, 1]
                fpr, tpr, _ = roc_curve(self.y_test, y_pred_proba)
                auc_score = auc(fpr, tpr)
                
                plt.plot(fpr, tpr, label=f'{model_name} (AUC = {auc_score:.3f})')
            
            plt.plot([0, 1], [0, 1], 'k--', label='Random Classifier')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('ROC Curves Comparison - RNA-only')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig('roc_curves_RNA.png', dpi=300, bbox_inches='tight')
            plt.show()
        
        # 4. Feature importance (for tree-based models) - MINT GREEN
        rf_model = self.results['Random Forest']['model']
        xgb_model = self.results['XGBoost']['model']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Random Forest importance
        rf_importance = pd.DataFrame({
            'feature': self.X_train_balanced.columns,
            'importance': rf_model.feature_importances_
        }).sort_values('importance', ascending=True).tail(20)
        
        sns.barplot(data=rf_importance, x='importance', y='feature', ax=ax1, color=MINT_GREEN)
        ax1.set_title('Random Forest - Top 20 Gene Importances (RNA-only)')
        
        # XGBoost importance
        xgb_importance = pd.DataFrame({
            'feature': self.X_train_balanced.columns,
            'importance': xgb_model.feature_importances_
        }).sort_values('importance', ascending=True).tail(20)
        
        sns.barplot(data=xgb_importance, x='importance', y='feature', ax=ax2, color=MINT_GREEN)
        ax2.set_title('XGBoost - Top 20 Gene Importances (RNA-only)')
        
        plt.tight_layout()
        plt.savefig('feature_importance_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        return True
    
    def detect_overfitting(self):
        """Analyze potential overfitting using cross-validation"""
        print("\n12. Analyzing overfitting with cross-validation...")
        
        cv_results = {}
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        for model_name, config in self.models.items():
            print(f"     Cross-validating {model_name}...")
            
            model = config['model']
            
            # Select appropriate data
            if config['requires_scaling']:
                X_data = self.X_train_scaled
            else:
                X_data = self.X_train_balanced
            
            train_scores = []
            val_scores = []
            
            for train_idx, val_idx in cv.split(X_data, self.y_train_balanced):
                X_train_fold = X_data.iloc[train_idx]
                X_val_fold = X_data.iloc[val_idx]
                y_train_fold = self.y_train_balanced[train_idx]
                y_val_fold = self.y_train_balanced[val_idx]
                
                # Train and evaluate
                model_copy = type(model)(**model.get_params())
                model_copy.fit(X_train_fold, y_train_fold)
                
                train_score = model_copy.score(X_train_fold, y_train_fold)
                val_score = model_copy.score(X_val_fold, y_val_fold)
                
                train_scores.append(train_score)
                val_scores.append(val_score)
            
            avg_train = np.mean(train_scores)
            avg_val = np.mean(val_scores)
            overfitting_gap = avg_train - avg_val
            
            cv_results[model_name] = {
                'Train_Accuracy': avg_train,
                'Validation_Accuracy': avg_val,
                'Overfitting_Gap': overfitting_gap,
                'Train_Std': np.std(train_scores),
                'Val_Std': np.std(val_scores)
            }
        
        cv_df = pd.DataFrame(cv_results).T.round(4)
        
        print("\n   Cross-Validation Results:")
        print(cv_df)
        
        # Visualize overfitting analysis
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Train vs Validation scores
        x_pos = np.arange(len(cv_df))
        width = 0.35
        
        ax1.bar(x_pos - width/2, cv_df['Train_Accuracy'], width, 
               label='Training Accuracy', alpha=0.8, color=MINT_GREEN_DARK)
        ax1.bar(x_pos + width/2, cv_df['Validation_Accuracy'], width, 
               label='Validation Accuracy', alpha=0.8, color=MINT_GREEN)
        
        ax1.set_xlabel('Models')
        ax1.set_ylabel('Accuracy')
        ax1.set_title('Training vs. Validation Accuracy - RNA-only')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(cv_df.index, rotation=45)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Overfitting gaps
        sns.barplot(x=cv_df.index, y=cv_df['Overfitting_Gap'], ax=ax2, color=MINT_GREEN)
        ax2.set_title('Overfitting Gap (Train - Validation) - RNA-only')
        ax2.set_ylabel('Accuracy Gap')
        ax2.axhline(y=0.05, color='red', linestyle='--', label='Warning Threshold (0.05)')
        ax2.axhline(y=0.1, color='red', linestyle='-', label='Critical Threshold (0.1)')
        ax2.legend()
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig('overfitting_analysis_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Overfitting warnings
        print("\n   Overfitting Analysis:")
        for model_name, results in cv_results.items():
            gap = results['Overfitting_Gap']
            if gap > 0.1:
                print(f"     WARNING: {model_name} shows significant overfitting (gap: {gap:.3f})")
            elif gap > 0.05:
                print(f"     CAUTION: {model_name} shows moderate overfitting (gap: {gap:.3f})")
            else:
                print(f"     GOOD: {model_name} shows minimal overfitting (gap: {gap:.3f})")
        
        self.cv_results = cv_df
        return cv_df
    
    def generate_final_report(self):
        """Generate comprehensive final report"""
        print("\n13. Generating final report...")
        
        # Create summary report
        report = {
            'Dataset_Info': {
                'Modality': 'scRNA-seq ONLY',
                'Total_Samples': len(self.y),
                'Total_Features_Original': self.X.shape[1] if hasattr(self, 'X') else 'N/A',
                'Total_Features_Final': self.X_engineered.shape[1],
                'Cell_Types': len(self.class_names),
                'Class_Names': list(self.class_names)
            },
            'Class_Distribution': self.imbalance_df.to_dict('records'),
            'Balancing_Applied': self.needs_balancing,
            'Model_Performance': self.metrics_df.to_dict('records'),
            'Cross_Validation': self.cv_results.to_dict() if hasattr(self, 'cv_results') else None
        }
        
        # Find best model
        best_model_idx = self.metrics_df['F1_Score'].idxmax()
        best_model_name = self.metrics_df.loc[best_model_idx, 'Model']
        best_f1_score = self.metrics_df.loc[best_model_idx, 'F1_Score']
        
        print(f"\n=== FINAL RESULTS (RNA-ONLY) ===")
        print(f"Best Model: {best_model_name}")
        print(f"Best F1-Score: {best_f1_score:.4f}")
        print(f"Dataset: {len(self.y)} cells, {self.X_engineered.shape[1]} features (genes)")
        print(f"Cell Types: {len(self.class_names)}")
        print(f"Class Balancing Applied: {self.needs_balancing}")
        
        # Save results
        with open('ml_pipeline_report_RNA.json', 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        # Save metrics tables
        self.metrics_df.to_csv('model_performance_summary_RNA.csv', index=False)
        
        # Save detailed results
        with pd.ExcelWriter('detailed_model_results_RNA.xlsx') as writer:
            self.metrics_df.to_excel(writer, sheet_name='Summary', index=False)
            if hasattr(self, 'cv_results'):
                self.cv_results.to_excel(writer, sheet_name='Cross_Validation')
            self.imbalance_df.to_excel(writer, sheet_name='Class_Distribution', index=False)
            
            for model_name, class_metrics in self.detailed_results.items():
                class_metrics.to_excel(writer, sheet_name=f'{model_name}_Details')
        
        print("\nFiles Generated:")
        print("  - ml_pipeline_report_RNA.json")
        print("  - model_performance_summary_RNA.csv")
        print("  - detailed_model_results_RNA.xlsx")
        print("  - Multiple visualization PNG files")
        
        return report
    
    def simple_feature_heatmap(self):
        """Simple feature importance heatmap using only tree-based models"""
        print("\n14. Creating simple feature importance heatmap...")
        
        # Only use models that definitely have feature_importances_
        tree_models = {}
        for model_name, results in self.results.items():
            if hasattr(results['model'], 'feature_importances_'):
                tree_models[model_name] = results['model'].feature_importances_
        
        if not tree_models:
            print("No tree-based models found for feature importance")
            return
        
        # Create simple heatmap
        importance_df = pd.DataFrame(tree_models, index=self.X_train_balanced.columns)
        top_features = importance_df.sum(axis=1).nlargest(20)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(importance_df.loc[top_features.index], annot=False, cmap='Greens')
        plt.title('Top 20 Gene Importance - Tree Models Only (RNA-only)')
        plt.tight_layout()
        plt.savefig('simple_feature_heatmap_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()

    def create_learning_curves(self):
        """Create learning curves to visualize overfitting"""
        print("\n15. Creating learning curves...")
        
        fig, axes = plt.subplots(1, len(self.models), figsize=(15, 5))
        if len(self.models) == 1:
            axes = [axes]
        
        for i, (model_name, config) in enumerate(self.models.items()):
            model = config['model']
            
            # Select appropriate data
            if config['requires_scaling']:
                X_data = self.X_train_scaled
            else:
                X_data = self.X_train_balanced
            
            # Generate learning curves
            train_sizes, train_scores, val_scores = learning_curve(
                model, X_data, self.y_train_balanced,
                cv=3, n_jobs=-1, train_sizes=np.linspace(0.1, 1.0, 8),
                scoring='accuracy'
            )
            
            train_mean = np.mean(train_scores, axis=1)
            train_std = np.std(train_scores, axis=1)
            val_mean = np.mean(val_scores, axis=1)
            val_std = np.std(val_scores, axis=1)
            
            axes[i].plot(train_sizes, train_mean, 'o-', color=MINT_GREEN_DARK, label='Training')
            axes[i].fill_between(train_sizes, train_mean - train_std, train_mean + train_std, alpha=0.1, color=MINT_GREEN_DARK)
            
            axes[i].plot(train_sizes, val_mean, 'o-', color=MINT_GREEN, label='Validation')
            axes[i].fill_between(train_sizes, val_mean - val_std, val_mean + val_std, alpha=0.1, color=MINT_GREEN)
            
            axes[i].set_title(f'{model_name} - RNA-only')
            axes[i].set_xlabel('Training Set Size')
            axes[i].set_ylabel('Accuracy')
            axes[i].legend()
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('learning_curves_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()

    def create_performance_radar(self):
        """Create radar chart comparing all models"""
        print("\n16. Creating performance radar chart...")
        
        # Prepare data
        metrics = ['Accuracy', 'Precision', 'Recall', 'F1_Score', 'AUC']
        models = self.metrics_df['Model'].tolist()
        
        # Number of variables
        N = len(metrics)
        
        # Angle for each axis
        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]  # Complete the circle
        
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        
        colors = [MINT_GREEN, MINT_GREEN_DARK, '#66CDAA']
        
        for i, model in enumerate(models):
            model_data = self.metrics_df[self.metrics_df['Model'] == model]
            values = [model_data[metric].iloc[0] for metric in metrics]
            values += values[:1]  # Complete the circle
            
            ax.plot(angles, values, 'o-', linewidth=2, label=model, color=colors[i % len(colors)])
            ax.fill(angles, values, alpha=0.25, color=colors[i % len(colors)])
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(metrics)
        ax.set_ylim(0, 1)
        ax.set_title('Model Performance Comparison (RNA-only)', size=16, pad=20)
        ax.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
        ax.grid(True)
        
        plt.savefig('performance_radar_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()

    def analyze_feature_distributions(self):
        """Analyze feature value distributions across cell types"""
        print("\n17. Creating feature distribution analysis...")
        
        # Get top features from Random Forest
        rf_model = self.results['Random Forest']['model']
        top_features_idx = np.argsort(rf_model.feature_importances_)[-6:]  # Top 6 features
        top_features = self.X_train_balanced.columns[top_features_idx]
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        for i, feature in enumerate(top_features):
            # Create DataFrame for plotting
            plot_data = pd.DataFrame({
                'Feature_Value': self.X_train_balanced[feature],
                'Cell_Type': [self.class_names[y] for y in self.y_train_balanced]
            })
            
            # Create violin plot - mint green palette
            sns.violinplot(data=plot_data, x='Cell_Type', y='Feature_Value', ax=axes[i], color=MINT_GREEN)
            axes[i].set_title(f'{feature} (RNA-only)')
            axes[i].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig('feature_distributions_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()

    def create_class_separation_plot(self):
        """Create PCA plot showing class separation"""
        print("\n18. Creating class separation visualization...")
        
        # Apply PCA
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(self.X_train_balanced)
        
        # Create DataFrame for plotting
        pca_df = pd.DataFrame({
            'PC1': X_pca[:, 0],
            'PC2': X_pca[:, 1],
            'Cell_Type': [self.class_names[y] for y in self.y_train_balanced]
        })
        
        plt.figure(figsize=(12, 8))
        
        # Create scatter plot with green color palette
        colors_map = plt.cm.Greens(np.linspace(0.4, 0.9, len(self.class_names)))
        
        for idx, cell_type in enumerate(self.class_names):
            mask = pca_df['Cell_Type'] == cell_type
            plt.scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], 
                       label=cell_type, alpha=0.6, s=50, color=colors_map[idx])
        
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        plt.title('Cell Type Separation in PCA Space (RNA-only)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('class_separation_pca_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()

    def create_gene_expression_heatmap(self):
        """Create comprehensive gene expression heatmap"""
        print("\n19. Creating gene expression heatmap...")
        
        # Get top variable genes
        gene_vars = self.X_train_balanced.var().sort_values(ascending=False)
        top_genes = gene_vars.head(50).index  # Top 50 most variable genes
        
        # Create expression matrix
        expr_matrix = self.X_train_balanced[top_genes].copy()
        
        # Add cell type labels
        cell_type_labels = [self.class_names[y] for y in self.y_train_balanced]
        expr_matrix['Cell_Type'] = cell_type_labels
        
        # Sort by cell type for better visualization
        expr_matrix = expr_matrix.sort_values('Cell_Type')
        
        # Separate labels and data
        plot_labels = expr_matrix['Cell_Type']
        plot_data = expr_matrix.drop('Cell_Type', axis=1)
        
        # Create color map for cell types
        unique_types = sorted(self.class_names)
        type_colors = dict(zip(unique_types, plt.cm.Set3(np.linspace(0, 1, len(unique_types)))))
        row_colors = [type_colors[ct] for ct in plot_labels]
        
        # Create clustermap
        plt.figure(figsize=(14, 10))
        g = sns.clustermap(plot_data.T, 
                          cmap='RdYlGn',
                          col_colors=row_colors,
                          figsize=(14, 10),
                          cbar_kws={'label': 'Expression Level'},
                          yticklabels=True,
                          xticklabels=False,
                          col_cluster=True,
                          row_cluster=True)
        
        g.fig.suptitle('Gene Expression Heatmap - Top 50 Variable Genes (RNA-only)', 
                      y=0.98, fontsize=14)
        
        # Add legend for cell types
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=type_colors[ct], label=ct) for ct in unique_types]
        g.ax_heatmap.legend(handles=legend_elements, 
                           bbox_to_anchor=(1.02, 0.5), 
                           loc='center left',
                           title='Cell Types')
        
        plt.savefig('gene_expression_heatmap_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"   Heatmap created with {len(top_genes)} most variable genes")

    def identify_marker_genes_per_celltype(self):
        """Identify and visualize marker genes for each cell type - MINT GREEN BARS"""
        print("\n20. Identifying marker genes per cell type...")
        
        # Calculate mean expression per cell type for all genes
        marker_genes = {}
        
        for class_idx, cell_type in enumerate(self.class_names):
            # Get cells of this type
            cell_mask = self.y_train_balanced == class_idx
            
            # Calculate mean expression in this cell type
            mean_in_type = self.X_train_balanced.loc[cell_mask].mean()
            
            # Calculate mean expression in other cell types
            mean_in_others = self.X_train_balanced.loc[~cell_mask].mean()
            
            # Calculate fold change
            fold_change = mean_in_type / (mean_in_others + 1e-10)
            
            # Get top 10 marker genes (highest fold change)
            top_markers = fold_change.nlargest(10)
            marker_genes[cell_type] = top_markers
        
        # Visualize as bar plots - MINT GREEN
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        for idx, (cell_type, markers) in enumerate(marker_genes.items()):
            if idx < len(axes):
                # Create bar plot with MINT GREEN
                markers.plot(kind='barh', ax=axes[idx], color=MINT_GREEN, edgecolor=MINT_GREEN_DARK, linewidth=1.5)
                axes[idx].set_title(f'Top 10 Marker Genes: {cell_type}', fontsize=12, fontweight='bold')
                axes[idx].set_xlabel('Fold Change (vs other cell types)', fontsize=10)
                axes[idx].set_ylabel('Gene', fontsize=10)
                axes[idx].grid(True, alpha=0.3, axis='x')
                
                # Add value labels
                for i, v in enumerate(markers.values):
                    axes[idx].text(v, i, f' {v:.2f}', va='center', fontsize=9)
        
        plt.suptitle('Per-Cell-Type Marker Genes (RNA-only)', fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig('per_celltype_marker_genes_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print marker genes
        print("\n   Top Marker Genes per Cell Type:")
        for cell_type, markers in marker_genes.items():
            print(f"\n   {cell_type}:")
            print(f"   {markers.head(5).to_dict()}")
        
        self.marker_genes = marker_genes
        return marker_genes

    def create_gene_expression_network(self):
        """Create gene-gene expression network with cell types (replaces TF network)"""
        print("\n21. Creating gene expression network...")
        
        # Get top features from Random Forest
        rf_model = self.results['Random Forest']['model']
        top_features_idx = np.argsort(rf_model.feature_importances_)[-15:]  # Top 15 genes
        top_features = self.X_train_balanced.columns[top_features_idx].tolist()
        
        # Calculate correlations between top genes
        feature_data = self.X_train_balanced[top_features]
        correlation_matrix = feature_data.corr()
        
        # Create network graph
        G = nx.Graph()
        
        # Add nodes (genes)
        for feature in top_features:
            G.add_node(feature, node_type='gene')
        
        # Add cell type nodes
        for cell_type in self.class_names:
            G.add_node(cell_type, node_type='cell_type')
        
        # Add edges between highly correlated genes (correlation > 0.5)
        for i, gene1 in enumerate(top_features):
            for j, gene2 in enumerate(top_features):
                if i < j:  # Avoid duplicate edges
                    corr = correlation_matrix.loc[gene1, gene2]
                    if abs(corr) > 0.5:  # Strong correlation
                        G.add_edge(gene1, gene2, weight=abs(corr), edge_type='gene-gene')
        
        # Add edges between genes and cell types based on mean expression
        for gene in top_features:
            for class_idx, cell_type in enumerate(self.class_names):
                # Calculate mean expression of this gene in this cell type
                cell_mask = self.y_train_balanced == class_idx
                if cell_mask.sum() > 0:
                    mean_expr = self.X_train_balanced.loc[cell_mask, gene].mean()
                    overall_mean = self.X_train_balanced[gene].mean()
                    
                    # If expression is significantly higher than average, add edge
                    if mean_expr > overall_mean + self.X_train_balanced[gene].std():
                        G.add_edge(gene, cell_type, weight=mean_expr, edge_type='gene-celltype')
        
        # Create visualization
        plt.figure(figsize=(14, 10))
        
        # Position nodes using spring layout
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # Draw nodes
        gene_nodes = [node for node, data in G.nodes(data=True) if data['node_type'] == 'gene']
        celltype_nodes = [node for node, data in G.nodes(data=True) if data['node_type'] == 'cell_type']
        
        nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color=MINT_GREEN, 
                              node_size=300, alpha=0.8, label='Genes')
        nx.draw_networkx_nodes(G, pos, nodelist=celltype_nodes, node_color='lightcoral', 
                              node_size=500, alpha=0.8, label='Cell Types')
        
        # Draw edges
        gene_edges = [(u, v) for u, v, data in G.edges(data=True) if data['edge_type'] == 'gene-gene']
        celltype_edges = [(u, v) for u, v, data in G.edges(data=True) if data['edge_type'] == 'gene-celltype']
        
        nx.draw_networkx_edges(G, pos, edgelist=gene_edges, edge_color='gray', 
                              alpha=0.5, width=1, style='dashed')
        nx.draw_networkx_edges(G, pos, edgelist=celltype_edges, edge_color=MINT_GREEN_DARK, 
                              alpha=0.7, width=2)
        
        # Draw labels
        labels = {}
        for node in G.nodes():
            if len(node) > 10:
                labels[node] = node[:10] + '...'
            else:
                labels[node] = node
        
        nx.draw_networkx_labels(G, pos, labels, font_size=8)
        
        plt.title('Gene Expression Network - RNA-only\n(Mint Green: Genes, Red: Cell Types, Dashed: Gene correlations, Solid: Gene-CellType associations)', 
                 fontsize=12)
        plt.legend()
        plt.axis('off')
        plt.tight_layout()
        plt.savefig('gene_expression_network_RNA.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print network statistics
        print(f"  Network created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        print(f"  Genes: {len(gene_nodes)}")
        print(f"  Cell types: {len(celltype_nodes)}")
        
        return G

    def run_complete_pipeline(self):
        """Execute the complete pipeline"""
        print("=== Starting Complete ML Pipeline (RNA-ONLY) ===")
        
        # Execute all steps
        steps = [
            self.load_data,
            self.analyze_class_imbalance,
            self.feature_engineering,
            self.prepare_train_test_split,
            self.handle_class_imbalance,
            self.scale_features,
            self.define_models,
            self.train_and_evaluate_models,
            self.calculate_metrics,
            self.create_detailed_evaluation_table,
            self.create_visualizations,
            self.detect_overfitting,
            self.generate_final_report,
            self.simple_feature_heatmap,
            self.create_learning_curves,
            self.create_performance_radar,
            self.analyze_feature_distributions,
            self.create_class_separation_plot,
            self.create_gene_expression_heatmap,
            self.identify_marker_genes_per_celltype,
            self.create_gene_expression_network,
        ]
        
        for i, step in enumerate(steps, 1):
            try:
                step()
                print(f"✓ Step {i}/{len(steps)} completed successfully")
            except Exception as e:
                print(f"✗ Step {i}/{len(steps)} failed: {e}")
                import traceback
                traceback.print_exc()
                return False
        
        print("\n=== Pipeline Completed Successfully! ===")
        return True

# ================================================================================================
# USAGE AND EXECUTION
# ================================================================================================

def main():
    """Main execution function"""
    print("="*80)
    print("Starting scRNA-seq Cell Type Classification Pipeline (RNA-ONLY)")
    print("="*80)
    
    # Initialize pipeline with correct path
    pipeline = scRNAMLPipeline(data_dir="/home/toheeb.jumah/Codeathon/python_ready_data")
    
    # Run complete pipeline
    success = pipeline.run_complete_pipeline()
    
    if success:
        print("\n" + "="*80)
        print("Pipeline completed successfully!")
        print("="*80)
        print("\n📊 Key Results:")
        print(f"- Best Model: {pipeline.metrics_df.loc[pipeline.metrics_df['F1_Score'].idxmax(), 'Model']}")
        print(f"- Best F1-Score: {pipeline.metrics_df['F1_Score'].max():.4f}")
        print(f"- Models Evaluated: {len(pipeline.models)}")
        print(f"- Features Used: {pipeline.X_engineered.shape[1]} (genes)")
        print(f"- Modality: scRNA-seq ONLY")
        
        # Show final metrics table
        print("\n📈 Final Model Comparison (RNA-only):")
        print(pipeline.metrics_df.round(4).to_string(index=False))
        
        print("\n📁 Generated Files:")
        print("   - All PNG visualizations with '_RNA' suffix")
        print("   - ml_pipeline_report_RNA.json")
        print("   - model_performance_summary_RNA.csv")
        print("   - detailed_model_results_RNA.xlsx")
        
    else:
        print("\n❌ Pipeline failed. Check error messages above!")

# Execute if run directly
if __name__ == "__main__":
    main()
