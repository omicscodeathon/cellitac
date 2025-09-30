
# Cell-Type Classification Analysis Summary Report
**Generated on:** 2025-09-29 11:11:20

## Dataset Summary
- **Total Cells:** 1,455
- **Total Features:** 7,000
- **Cell Types:** 4
- **Training Cells:** 1,017
- **Test Cells:** 438

## Model Performance Summary
              Model  Train_Accuracy  Test_Accuracy  CV_Mean   CV_Std  Overfitting_Score   Overfitting_Status
      Random Forest        0.998033       0.947489 0.951835 0.012535           0.050545 Moderate Overfitting
  Gradient Boosting        1.000000       0.947489 0.947880 0.010621           0.052511 Moderate Overfitting
Logistic Regression        1.000000       0.974886 0.971487 0.005727           0.025114             Good Fit
                SVM        1.000000       0.917808 0.926263 0.005248           0.082192 Moderate Overfitting
     Neural Network        1.000000       0.968037 0.963629 0.005834           0.031963             Good Fit

## Best Performing Models
**Highest Test Accuracy:** Logistic Regression (0.9749)
**Most Stable (Lowest CV Std):** SVM (±0.0052)
**Best Fit (Lowest Overfitting):** Logistic Regression (0.0251)

## Timing Analysis Summary
              Model  Training_Time_s  Prediction_Time_s  Test_Accuracy  Efficiency_Score  Speed_Rank  Accuracy_Rank  Efficiency_Rank
      Random Forest         0.344961           0.109734       0.947489          2.746653           1              3                1
  Gradient Boosting       388.737365           0.109868       0.947489          0.002437           5              3                5
Logistic Regression         1.387288           0.077145       0.974886          0.702728           2              1                2
                SVM        15.067777           4.751080       0.917808          0.060912           4              5                4
     Neural Network         2.810880           0.106075       0.968037          0.344389           3              2                3

**Fastest Training:** Random Forest (0.34s)
**Most Efficient:** Random Forest (2.747)

## Cell Type Analysis
**Cell Types Analyzed:** Monocyte, T_cells, NK_cell, B_cell

## Files Generated
### Models
- Trained models: `results/models/`
- Label encoder: `results/models/label_encoder.pkl`
- Feature scaler: `results/models/feature_scaler.pkl`

### Analysis Results
- Model performance: `results/analysis/model_performance_summary.csv`
- Feature importance: `results/analysis/*_feature_importance.csv`
- Cell-type features: `results/analysis/celltype_*_features.csv`
- TF activity: `results/analysis/tf_activity_analysis.csv`
- Top TFs: `results/analysis/top_tfs_*.csv`
- Timing analysis: `results/analysis/model_timing_analysis.csv`

### Visualizations
- Performance plots: `results/plots/model_performance_comparison.png`
- Feature heatmaps: `results/plots/*_heatmap.png`
- TF network: `results/plots/tf_celltype_network.png`
- TF activity: `results/plots/tf_activity_dotplot.png`

## Recommendations
1. **Best Overall Model:** Logistic Regression for highest accuracy
2. **Production Model:** Random Forest for best efficiency
3. **Overfitting Concerns:** Monitor 3 models showing overfitting

## Next Steps
1. Validate findings with additional datasets
2. Investigate top TF features for biological relevance
3. Consider ensemble methods for improved performance
4. Optimize hyperparameters for best-performing models
