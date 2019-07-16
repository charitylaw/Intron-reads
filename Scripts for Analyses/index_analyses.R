library(index)
library(magrittr)
library(stringr)

## Helper functions ----
plot_index_barplot <- function(index_output) {
    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 7:8, 5:6, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    categories <- index_output$category
    categories <- factor(
        categories,
        levels = c("Mixed+-", "Mixed-+", "Intron-", "Intron+", "Exon-", "Exon+", "-", "+", "")
    )

    category_counts <- sapply(
        setdiff(levels(categories), ""),
        function(x) {
            sum(categories == x)
        }
    )
    count_o <- order(category_counts)
    category_counts <- category_counts[count_o]
    bar_col <- bar_col[count_o]

    barplot(
        category_counts,
        col = bar_col,
        cex.axis = 1,
        cex.lab = 1,
        las = 2,
        main = "index Categories"
    )
}

plot_index_tstat <- function(index_output) {
    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 7:8, 5:6, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    categories <- index_output$category
    categories <- factor(
        categories,
        levels = c("Mixed+-", "Mixed-+", "Intron-", "Intron+", "Exon-", "Exon+", "-", "+", "")
    )
    point_col <- bar_col[categories]

    category_counts <- sapply(
        setdiff(levels(categories), ""),
        function(x) {
            sum(categories == x)
        }
    )
    count_o <- order(category_counts)
    category_counts <- category_counts[count_o]
    bar_col <- bar_col[count_o]

    plot(
        index_output$tops$intron$t,
        index_output$tops$exon$t,
        pch = 20,
        cex.axis = 1,
        cex.lab = 1,
        col = point_col,
        main = "t-statistic",
        xlab = "Intron t-statistic",
        ylab = "Exon t-statistic"
    )
    abline(a = 0, b = 1, col = "#BBBBBB", lwd = 2)
}

plot_index_logfc <- function(index_output) {
    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 7:8, 5:6, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    categories <- index_output$category
    categories <- factor(
        categories,
        levels = c("Mixed+-", "Mixed-+", "Intron-", "Intron+", "Exon-", "Exon+", "-", "+", "")
    )
    point_col <- bar_col[categories]

    category_counts <- sapply(
        setdiff(levels(categories), ""),
        function(x) {
            sum(categories == x)
        }
    )
    count_o <- order(category_counts)
    category_counts <- category_counts[count_o]
    bar_col <- bar_col[count_o]

    plot(
        index_output$tops$intron$logFC,
        index_output$tops$exon$logFC,
        pch = 20,
        cex.axis = 1,
        cex.lab = 1,
        col = point_col,
        main = "logFC",
        xlab = "Intron logFC",
        ylab = "Exon logFC"
    )
    abline(a = 0, b = 1, col = "#BBBBBB", lwd = 2)
}

plot_index_boxplot <- function(index_output) {
    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 7:8, 5:6, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    categories <- index_output$category
    categories <- factor(
        categories,
        levels = c("Mixed+-", "Mixed-+", "Intron-", "Intron+", "Exon-", "Exon+", "-", "+", "")
    )

    category_counts <- sapply(
        setdiff(levels(categories), ""),
        function(x) {
            sum(categories == x)
        }
    )
    count_o <- order(category_counts)
    category_counts <- category_counts[count_o]
    bar_col <- bar_col[count_o]

    merged_categories <- categories
    levels(merged_categories) <- c("mixed", "mixed", "intron", "intron", "exon", "exon", "both", "both", "0")
    boxplot_cols <- RColorBrewer::brewer.pal(9, "Paired")[c(3:4, 7:8, 5:6, 1:2, 9)]
    boxplot_cols <- boxplot_cols[c(2,4,6,8,9)]

    in_len_split <- split(sqrt(index_output$dges$intron$genes$Length), merged_categories)
    in_len_mean <- sapply(in_len_split, median)

    in_len_o <- order(in_len_mean)
    in_len_split <- in_len_split[in_len_o]
    boxplot_cols <- boxplot_cols[in_len_o]

    boxplot(
        in_len_split,
        outline = FALSE,
        border = boxplot_cols,
        cex.axis = 1,
        cex.lab = 1,
        lwd = 2,
        las = 2,
        yaxt = "n",
        main = "Total intron length (in thousands)"
    )

    ax <- ((0:7)*100)
    ax[2] <- sqrt(20000)
    axis(side = 2, at = ax, labels = (ax^2)/1e3, las = 2)
    axis(side = 2, at = sqrt(median(index_output$dges$exon$genes$Length)), labels = "exon \n median", las = 2, font = 3, col.axis = "black", col.ticks = "black")
}

## Taken from https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region = "figure", pos = "topleft", cex = NULL, ...) {
    region <- match.arg(region, c("figure", "plot", "device"))
    pos <- match.arg(
        pos,
        c(
            "topleft",
            "top",
            "topright",
            "left",
            "center",
            "right",
            "bottomleft",
            "bottom",
            "bottomright"
        )
    )

    if (region %in% c("figure", "device")) {
        ds <- dev.size("in")
        # xy coordinates of device corners in user coordinates
        x <- grconvertX(c(0, ds[1]), from = "in", to = "user")
        y <- grconvertY(c(0, ds[2]), from = "in", to = "user")

        # fragment of the device we use to plot
        if (region == "figure") {
            # account for the fragment of the device that
            # the figure is using
            fig <- par("fig")
            dx <- (x[2] - x[1])
            dy <- (y[2] - y[1])
            x <- x[1] + dx * fig[1:2]
            y <- y[1] + dy * fig[3:4]
        }
    }

    # much simpler if in plotting region
    if (region == "plot") {
        u <- par("usr")
        x <- u[1:2]
        y <- u[3:4]
    }

    sw <- strwidth(text, cex = cex)
    sh <- strheight(text, cex = cex)

    x1 <- switch(
        pos,
        topleft     = x[1] + sw,
        left        = x[1] + sw,
        bottomleft  = x[1] + sw,
        top         = (x[1] + x[2]) / 2,
        center      = (x[1] + x[2]) / 2,
        bottom      = (x[1] + x[2]) / 2,
        topright    = x[2] - sw,
        right       = x[2] - sw,
        bottomright = x[2] - sw
    )

    y1 <- switch(
        pos,
        topleft     = y[2] - sh,
        top         = y[2] - sh,
        topright    = y[2] - sh,
        left        = (y[1] + y[2]) / 2,
        center      = (y[1] + y[2]) / 2,
        right       = (y[1] + y[2]) / 2,
        bottomleft  = y[1] + sh,
        bottom      = y[1] + sh,
        bottomright = y[1] + sh
    )

    old.par <- par(xpd = NA)
    on.exit(par(old.par))

    text(x1, y1, text, cex = cex, ...)
    return(invisible(c(x, y)))
}

# Main computation ----
index_immune <- local({
    exon <- readRDS("Immune/ImmuneExons.Rds")
    intron <- readRDS("Immune/ImmuneIntrons.Rds")
    genebody <- readRDS("Immune/ImmuneGenebody.Rds")
    group <- readRDS("Immune/ImmuneGroups.Rds")

    intron$genes$Length <- genebody$genes$Length - exon$genes$Length + 1

    design <- model.matrix(~ 0 + group) %>%
        set_colnames(colnames(.) %>% str_remove("group"))

    get_pairwise_contrasts <- function(design, cols = ncol(design)) {
        group_names <- colnames(design)

        contrasts <- list()
        contrast_names <- character()
        for (i in 1:(cols - 1)) {
            for (j in (i + 1):cols) {
                contr <- numeric(length = cols)
                contr[i] <- 1
                contr[j] <- -1
                contrast_name <- glue::glue("{group_names[i]} vs {group_names[j]}")
                contrasts <- append(contrasts, list(contr))
                contrast_names <- append(contrast_names, contrast_name)
            }
        }
        setNames(contrasts, contrast_names)
    }

    pairwise_contrasts <- get_pairwise_contrasts(design)

    index_analysis(exon, intron, group = group, design = design, contrast = pairwise_contrasts[["Monocytes vs Neutrophils"]], p.value = 0.01)
})

index_cell_line <- local({
    exon <- readRDS("CellLine/CellLineExons.Rds")
    intron <- readRDS("CellLine/CellLineIntrons.Rds")
    group <- readRDS("CellLine/CellLineGroups.Rds")

    index_analysis(exon, intron, group, p.value = 0.01)
})

# Output plotting ----
local({
    pdf("index_figure_4.pdf", height = 6, width = 11)
    par(mfrow = c(2, 4))

    plot_index_barplot(index_cell_line)
    fig_label("(a)", cex = 1.5)
    plot_index_tstat(index_cell_line)
    fig_label("(b)", cex = 1.5)
    plot_index_logfc(index_cell_line)
    fig_label("(c)", cex = 1.5)
    plot_index_boxplot(index_cell_line)
    fig_label("(d)", cex = 1.5)
    plot_index_barplot(index_immune)
    fig_label("(e)", cex = 1.5)
    plot_index_tstat(index_immune)
    fig_label("(f)", cex = 1.5)
    plot_index_logfc(index_immune)
    fig_label("(g)", cex = 1.5)
    plot_index_boxplot(index_immune)
    fig_label("(h)", cex = 1.5)

    par(mfrow = c(1, 1))
    dev.off()
})
