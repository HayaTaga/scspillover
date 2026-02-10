dir.create('data/nonproprietary', recursive = TRUE, showWarnings = FALSE)

data(california_smoking, package = 'scspill')
write.csv(california_smoking$panel, 'data/nonproprietary/california_panel.csv', row.names = FALSE)
write.csv(california_smoking$w, 'data/nonproprietary/california_w_vector.csv', row.names = FALSE)
write.csv(california_smoking$W, 'data/nonproprietary/california_W_matrix.csv', row.names = FALSE)

load('data/sudan_secession.rda')
write.csv(sudan_secession$panel_df, 'data/nonproprietary/sudan_panel.csv', row.names = FALSE)
write.csv(sudan_secession$w, 'data/nonproprietary/sudan_w_vector.csv', row.names = FALSE)
write.csv(sudan_secession$W, 'data/nonproprietary/sudan_W_matrix.csv', row.names = FALSE)

message('Exported CSV copies into data/nonproprietary/.')
