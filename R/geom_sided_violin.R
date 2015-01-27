library(proto)

GeomSidedViolin <- proto(ggplot2:::GeomViolin, {
    draw <- function(., data, ...) {
        # Find the points for the line to go all the way around
        data <- transform(data, xminv = x - violinwidth * (side == 0) * (x-xmin),
                          xmaxv = x + violinwidth * (side == 1) * (xmax-x))

        newdata <- rbind(arrange(transform(data, x = xminv), y),
                         arrange(transform(data, x = xmaxv), -y))

        newdata <- rbind(newdata, newdata[1,])

        ggname(.$my_name(), GeomPolygon$draw(newdata, ...))
    }
})


#' a violin plot with two sides
#'
#' A sided violin plot is good for comparing two distributions across
#' many sets, such as comparing microarray and RNA-Seq for each set
#'
#' @export
geom_sided_violin <- function (mapping = NULL, data = NULL, stat = "ydensity",
                              position = "dodge", trim = TRUE,
                              scale = "area", ...)
{
    GeomSidedViolin$new(mapping = mapping, data = data, stat = stat,
                        position = position, trim = trim, scale = scale, ...)
}
