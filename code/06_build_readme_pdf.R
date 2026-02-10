input_file <- "README_submit.md"
output_file <- "README.pdf"

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

lines <- readLines(input_file, warn = FALSE)
lines_per_page <- 55L
line_height <- seq(0.98, 0.02, length.out = lines_per_page)

pdf(output_file, width = 8.5, height = 11)
for (i in seq(1L, length(lines), by = lines_per_page)) {
  page_lines <- lines[i:min(i + lines_per_page - 1L, length(lines))]
  plot.new()
  par(mar = c(0, 0, 0, 0))
  text(
    x = rep(0.02, length(page_lines)),
    y = line_height[seq_along(page_lines)],
    labels = page_lines,
    adj = c(0, 1),
    cex = 0.7,
    family = "mono"
  )
}
dev.off()

message(sprintf("Wrote %s", output_file))
