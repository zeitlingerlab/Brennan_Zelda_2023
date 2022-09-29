library(lattice)

normalize_matrix <- function(m, value) {
  m <- pmin(m / value, 1)
  m
}

draw_standard_heatmap <- function(m, title, center_vline=FALSE, normalize=TRUE, colors=colorRampPalette(c("#F6F6F6", "blue"))(32)) {
  
  if(normalize) {
    max.value <- quantile(m, 0.98, na.rm=TRUE)
    m <- normalize_matrix(m, max.value)
  }

  image(t(m[nrow(m):1, ]), col=colors, useRaster=TRUE, main=title, yaxt="n", xaxt="n")

  if(center_vline) {
    vline <- matrix(NA, nrow=nrow(m), ncol=ncol(m))
    
    m_center_col <- ceiling(ncol(m) / 2)
    m_line_width <- ceiling(ncol(m) * 0.01)
    m_center_range <- (m_center_col - m_line_width) : (m_center_col + m_line_width)
    vline[, m_center_range] <- matrix(1, nrow=nrow(vline), ncol=length(m_center_range))
    vline.colors <- paste0(colorRampPalette("gray")(1), "77")
    image(t(vline[nrow(vline):1,]), col=vline.colors, useRaster=TRUE, add=TRUE, yaxt="n", xaxt="n")
  }
}

rbind_matrices_with_blanks <- function(mlist, blank_row_count) {
  blanks <- matrix(NA, nrow=blank_row_count, ncol=ncol(mlist[[1]]))
  first_one <- TRUE
  m <- NULL
  for(i in 1:length(mlist)) {
    if(first_one == TRUE) {
      m <- mlist[[i]]
      first_one <- FALSE
    } else {
      m <- rbind(m, blanks, mlist[[i]])
    }
  }
  m
}

