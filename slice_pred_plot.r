slice_pred_plot <- function(data, df_estimate, slice, method) {
  library(ggplot2)
  if (slice == "NGP") {
    newdata <- data.frame(Cultivar = data$Cultivar,
                          NGP = data$NGP, 
                          MHG = median(data$MHG),
                          GY = data$GY)
     }else {newdata <- data.frame(Cultivar = data$Cultivar,
                          NGP = median(data$NGP),
                          MHG = data$MHG,
                          GY = data$GY)}
  
  predictions <- data.frame()
  # Loop over each row in 'newdata'
  for (i in 1:nrow(newdata)) {
    beta_df <- df_estimate[df_estimate$cultivar == newdata$Cultivar[i], ]
    pred <- beta_df$beta_0 + 
            beta_df$beta_NGP * newdata$NGP[i] + 
            beta_df$beta_MHG * newdata$MHG[i] + 
            beta_df$beta_NGP_MHG * newdata$NGP[i] * newdata$MHG[i] + 
            beta_df$beta_MHG2 * newdata$MHG[i]^2
    
    # Add the predictions to the 'predictions_bayes' data frame
    predictions <- rbind(predictions, data.frame(cultivar = newdata$Cultivar[i],
                                                  NGP = newdata$NGP[i],
                                                  MHG = newdata$MHG[i],
                                                  prediction = pred))
  }
  # Plot predictions against slice
  if (slice == "NGP"){
    ggplot() +
    geom_line(data = predictions, aes(x = NGP, y = prediction, color = cultivar)) +
    labs(title = paste(method, " Predictions (NGP slice)"),
         x = slice,
         y = "GY") +
    ylim(-1, 3.5) +
    xlim(0, 1) +
    theme_minimal() +
    theme(legend.position = "none")
  }
  else{
    ggplot() +
      geom_line(data = predictions, aes(x = MHG, y = prediction, color = cultivar)) +
      labs(title = paste(method, " Predictions (MHG slice)"),
        x = slice,
        y = "GY") +
      theme_minimal() +
      theme(legend.position = "none")
  }
}