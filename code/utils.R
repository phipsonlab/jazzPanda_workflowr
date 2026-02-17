# define figure and data path
library(here)
library(ggplot2)
library(patchwork)
library(dotenv)
# ------------------------------------------------------------------------------
# Define directories for raw and processed data used in the jazzPanda paper
# ------------------------------------------------------------------------------
dotenv::load_dot_env(file = here(".env"))
raw_dir <- Sys.getenv("RAW_DATA_DIR")

# Paths to specific datasets used in the analyses
rdata <- list(
    
    # ---------------- Xenium mouse brain ----------------
    # Replicate 1: FF multi-section dataset
    xmb_r1 = file.path(raw_dir, "Xenium_mouse_brain", 
                       "Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs/"),
    
    # Replicate 2: second multi-section dataset
    xmb_r2 = file.path(raw_dir, "Xenium_mouse_brain", 
                       "Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs/"),
    
    # Replicate 3: third multi-section dataset
    xmb_r3 = file.path(raw_dir, "Xenium_mouse_brain", 
                       "Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/"),
    
    # ---------------- Xenium human breast ----------------
    # Human breast cancer sample 1 (Xenium standard panel)
    xhb_r1 = file.path(raw_dir, "Xenium_human_breast_samples", 
                       "Xenium_hbreast_sample1/"),
    
    # Human breast cancer sample 2 (Xenium standard panel)
    xhb_r2 = file.path(raw_dir, "Xenium_human_breast_samples", 
                       "Xenium_hbreast_sample2/"),
    
    # ---------------- Xenium human lung ----------------
    # Human lung cancer dataset (single sample)
    xhl = file.path(raw_dir, "Xenium_human_lung_cancer/"),
    
    # ---------------- MERSCOPE human breast ----------------
    # Human breast cancer patient sample (MERSCOPE platform, external directory)
    mhb = file.path(raw_dir, "Merscope_human_breast_sample/"),
    
    # ---------------- CosMx human liver ----------------
    # Healthy and diseased liver samples (CosMx platform)
    cosmx_data = file.path(raw_dir, "CosMx_normal_and_diseased_liver_samples")

)

# ---------------- Processed / generated datasets ----------------
# Pre-computed / processed data (e.g. complexity benchmarks)
# from Zenodo/GitHub associated with the jazzPanda paper
gdata = Sys.getenv("PROCESSED_DATA_DIR")
# cluster colors
my_colors <- c(
    "#800000", "#f032e6","#ffe119", "#e6beff",
    "#e6194b", "#008080", "#0082c8", "#3cb44b",
    "#aa6e28", "#fabebe", "#d2f53c", "#911eb4",
    "#fffac8", "#f58231", "#46f0f0","#808080"
)

#############################################################################
# theme setting for figures
defined_theme <- theme(strip.text = element_text(size = rel(2)),
                       axis.line=element_blank(),
                       axis.ticks=element_blank(),
                       axis.text=element_blank(),
                       legend.position = "none",
                       axis.title=element_blank(),
                       panel.background=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.background=element_blank(),
                       panel.border = element_rect(fill=NA, colour = "black"))


#############################################################################
# Function to automatically determine hex bin size based on data density
# Function to automatically determine hex bin size based on data density
auto_hex_bin <- function(n, target_points_per_bin = 0.5, min_bins = 15, 
                         max_bins = 200) {
    if (n < 50) {
        return(min_bins)
    }
    k <- n / target_points_per_bin
    bins <- round(sqrt(k))
    
    # For highly expressed genes, allow more bins to resolve spatial structure
    # Scale up target more aggressively past a threshold
    if (n > 5000) {
        bins <- round(sqrt(n / 3))
    }
    
    return(min(max_bins, max(min_bins, bins)))
}


#############################################################################
# Function to calculate the cumulative average correlation 
get_cmr_ma<- function(genes, cor_M, cl){
    curr_corrs = cor_M[genes, cl]
    mv_avgs = cumsum(curr_corrs)/(1:length(genes))
    return(mv_avgs)
}

#############################################################################
# A help function to load Xenium data
get_xenium_data<-function(path,mtx_name, trans_name="transcript_info.csv.gz",
                          cells_name="cell_info.csv.gz"){
    
    transcript_info <- as.data.frame(fread(paste(path, trans_name,sep="")))
    cell_info <- as.data.frame(fread(paste(path,cells_name,sep="")))
    
    data <- Read10X(data.dir = paste(path,mtx_name, sep=""))
    
    cm <- as.matrix(data$`Gene Expression`)
    
    r_codeword <- as.matrix(data$`Negative Control Codeword`)
    
    r_probe <- as.matrix(data$`Negative Control Probe`)
    # merge negative control genes and real genes
    cm_neg <- as.data.frame(rbind(r_probe, r_codeword))
    zero_cells <- colnames(cm)[colSums(cm)==0]
    
    cm = cm[, setdiff(colnames(cm),zero_cells)]
    transcript_info$x <- as.numeric(transcript_info$x_location)
    transcript_info$y <- as.numeric(transcript_info$y_location)
    
    cell_info$x <- as.numeric(cell_info$x_centroid)
    cell_info$y <- as.numeric(cell_info$y_centroid)
    
    return (list(cm = cm, cm_neg=cm_neg, zero_cells = zero_cells,
                 trans_info=transcript_info, cell_info=cell_info,
                 probe = row.names(r_probe),
                 codeword=row.names(r_codeword)))
    
}


#############################################################################
# Function to create one formatted summary line
make_sc_summary <- function(x, label) {
    s <- summary(x)
    sprintf("%s\n  Min: %.4f |  1Q: %.4f |  Median: %.4f |  Mean: %.4f |  3Q: %.4f |  Max: %.4f\n",
            label, s[1], s[2], s[3], s[4], s[5], s[6])
}


#############################################################################
# help function to plot UMAP
plot_umap_seu <- function(cluster_info, seu, file_prefix,
                          ct_nm,
                          my_colors,
                          out_dir = overview_PA, fig_w=1000) {
    
    
    # Build output file path
    out_file <- here(out_dir, paste0(file_prefix, "_UMAP", ".jpg"))
    # Save to PDF
    jpeg(out_file, width = fig_w, height = 1000, res=200)
    print(DimPlot(seu, reduction = "umap",split.by = "sample")+  ggtitle("") +
              labs(title = " ", x = "UMAP1", y = "UMAP2", fill=" ") +
              scale_color_manual(values = my_colors)+
              theme(legend.box.margin = margin(0, 0, 0, 2), 
                    plot.margin = margin(5, 5, 5, 5),     
                    legend.position = "right",
                    panel.border = element_blank(), 
                    panel.grid = element_blank(),
                    legend.text = element_text(size = 12),
                    panel.background = element_blank(),
                    axis.title = element_text(size=12)
              ) )
    dev.off()
}

#############################################################################
plot_top3_genes <- function(vis_df, fig_ratio = 1) {
    inters = unique(as.character(vis_df$feature_name))
    genes_plt <- ggplot(data = vis_df, aes(x = x, y = y)) +
        geom_hex(bins = auto_hex_bin(max(table(vis_df$feature_name)))) +
        facet_wrap(~text_label, nrow=1) +
        scale_fill_gradient(low="grey90", high="maroon4") + 
        # scale_fill_viridis_c(option = "turbo")+
        guides(fill = guide_colorbar(barheight = unit(0.2, "npc"), 
                                     barwidth  = unit(0.02, "npc")))+
        defined_theme +
        theme(
            legend.position = "right",
            aspect.ratio = fig_ratio,
            strip.background = element_rect(fill = NA, colour = NA),
            strip.text = element_text(size = rel(1.3)),
            strip.text.y.right = element_blank(),
            legend.key.width = unit(2, "cm"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.key.height = unit(2, "cm"),
            plot.margin = margin(0, 0, 0, 0)
        )
    
    if (length(inters) == 2) {
        genes_plt <- (genes_plt | plot_spacer()) +
            plot_layout(widths = c(2, 1), heights = c(1))
    } else if (length(inters) == 1) {
        genes_plt <- (genes_plt | plot_spacer() | plot_spacer()) +
            plot_layout(widths = c(1, 1, 1), heights = c(1))
    }
    
    return(genes_plt)
}