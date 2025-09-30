
# Cell-Type Classification Analysis Summary Report
**Generated on:** 2025-09-28 15:00:16

## Dataset Summary
- **Total Cells:** 1,463
- **Total Features:** 7,000
- **Cell Types:** 6
- **Training Cells:** 1,022
- **Test Cells:** 441

## Model Performance Summary
              Model  Train_Accuracy  Test_Accuracy  CV_Mean   CV_Std  Overfitting_Score   Overfitting_Status
      Random Forest        0.997065       0.938776 0.947154 0.012582           0.058289 Moderate Overfitting
  Gradient Boosting        1.000000       0.936508 0.941272 0.019371           0.063492 Moderate Overfitting
Logistic Regression        1.000000       0.965986 0.965744 0.008237           0.034014             Good Fit
                SVM        0.998043       0.911565 0.919742 0.012757           0.086478 Moderate Overfitting
     Neural Network        1.000000       0.965986 0.965744 0.008237           0.034014             Good Fit

## Best Performing Models
**Highest Test Accuracy:** Logistic Regression (0.9660)
**Most Stable (Lowest CV Std):** Logistic Regression (±0.0082)
**Best Fit (Lowest Overfitting):** Logistic Regression (0.0340)

## Timing Analysis Summary
              Model  Training_Time_s  Prediction_Time_s  Test_Accuracy  Efficiency_Score  Speed_Rank  Accuracy_Rank  Efficiency_Rank
      Random Forest         0.341454           0.117096       0.938776          2.749351           1              3                1
  Gradient Boosting       482.161740           0.113548       0.936508          0.001942           5              4                5
Logistic Regression         0.385279           0.079137       0.965986          2.507241           2              1                2
                SVM        17.446561           4.500669       0.911565          0.052249           4              5                4
     Neural Network         2.851532           0.108628       0.965986          0.338760           3              1                3

**Fastest Training:** Random Forest (0.34s)
**Most Efficient:** Random Forest (2.749)

## Cell Type Analysis
**Cell Types Analyzed:** T_cells, Monocyte, NK_cell, B_cell, HSC_-G-CSF, Pre-B_cell_CD34-

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
