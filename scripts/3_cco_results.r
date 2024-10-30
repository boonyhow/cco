# IMPORTS
libraries <- 'argparse,data.table,ggplot2,ggridges,ggrepel,glue,readxl,writexl,yaml,utils,dplyr,edgeR,tidyr,topGO,org.Hs.eg.db,tools'

# Split the string into a list of individual library names
libraries <- strsplit(libraries, ",")[[1]]

# Load each library using a loop, suppressing startup messages
for (name in libraries) {
  suppressPackageStartupMessages(library(name, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
}

set.seed(1)

# Define the class GeneExpressionAnalysis
GeneExpressionAnalysis <- setRefClass(
  "GeneExpressionAnalysis",
  
  fields = list(
    config = "list",
    marker_info = "data.table",
    rnaseq = "data.table",
    rnaseq_meta = "data.table",
    melted_rna_marker = "data.table",
    merged_rna_marker = "data.table",
    counts = "data.table",
    is_wide = "logical"  # NEW: Flag to indicate whether data is wide
  ),
  
  methods = list(
    
    initialize = function(config_file) {
      print("Loading configuration from YAML file...")
      config <<- read_yaml(config_file)
      print("Configuration loaded successfully.")
    },
    
    load_data = function() {
      print("Loading marker information...")
      marker_info <<- fread(glue("{config$paths$marker_merged_file}"))
      print(paste("Loaded", nrow(marker_info), "rows of marker info data."))
      
      marker_info <<- na.omit(marker_info)
      marker_info[, "V1" := NULL]


      # Check if the file is CSV or XLSX and load accordingly
      file_extension <- tools::file_ext(config$paths$rnaseq_file)
      print(file_extension)
      if (file_extension == "txt") {
        print("Loading RNA-seq data from TSV...")
        rnaseq <<- data.table(read.csv(config$paths$rnaseq_file, sep = '\t'))
        rnaseq_meta <<- rnaseq %>%
        dplyr::select(sample, day) %>%
        distinct()  # Retain unique rows for metadata

      # Convert the day column to a numeric day by extracting the day number
        rnaseq_meta <<- rnaseq_meta %>%
          mutate(day_num = as.numeric(gsub("day", "", day)))  # Strip "day" prefix and convert to numeric

      } else if (file_extension == "csv") {
        print("Loading RNA-seq data from CSV...")
        rnaseq <<- data.table(read.csv(config$paths$rnaseq_file, sep = ','))
        rnaseq_meta <<- rnaseq %>%
        dplyr::select(sample, day) %>%
        distinct()  # Retain unique rows for metadata

      # Convert the day column to a numeric day by extracting the day number
        rnaseq_meta <<- rnaseq_meta %>%
          mutate(day_num = as.numeric(gsub("day", "", day)))  # Strip "day" prefix and convert to numeric

      } else if (file_extension == "xlsx") {
        print("Loading RNA-seq data from XLSX...")
        rnaseq <<- data.table(read_excel(config$paths$rnaseq_file, sheet='log2TPM'))
         rnaseq_meta <<-  data.table(read_excel(config$paths$rnaseq_file, sheet='Samples'))
      } else {
        stop("Unsupported file type. Please provide a .txt, .csv or .xlsx file for RNA-seq data.")
      }

      rnaseq <<- rnaseq[gene_type == "protein_coding"]
      print(paste("Loaded", nrow(rnaseq), "rows of RNA-seq data."))

      print("Loading RNA-seq metadata...")

      print(paste("Loaded", nrow(rnaseq_meta), "rows of RNA-seq metadata."))
    },
    
    countTrues = function(vec) {
      vec <- as.numeric(vec)
      return(sum(vec))
    },
    
    calculate_cco_levels = function() {
      print("Calculating the cco_levels of true values in the marker info data...")
      marker_info$cco_levels <<- apply(marker_info[, 2:length(colnames(marker_info))], 1, .self$countTrues)
      print("cco_levels calculation completed.")
    },
    
    melt_rna = function() {
      print("Checking if RNA-seq data is in wide or long format...")

      # Check if 'SRR' or 'ENCL' exists in the column names
      if (any(grepl("*ENCL", colnames(rnaseq)))) {
          pattern <- "*ENCL"
          print(glue("Using '{pattern}' pattern for melting."))
          is_wide <<- TRUE
      } else if (any(grepl("*SRR", colnames(rnaseq)))) {
          pattern <- "*SRR"
          print(glue("Using '{pattern}' pattern for melting."))
          is_wide <<- TRUE
      } else {
          # If neither pattern exists, assume the table is already in long format
          print("No 'SRR' or 'ENCL' patterns found. Assuming the table is already in long format.")
          melted_rna_marker <<- rnaseq  # Assign the existing table as the long format
          melted_rna_marker$gene_id <- gsub("\\..*", "", melted_rna_marker$gene_id)
          print(paste("Assumed long RNA-seq data with", nrow(melted_rna_marker), "rows."))
          is_wide <<- FALSE  # Data is long
          return()  # Exit the function if data is already long
      }

      # If data is wide, melt it
      print("Melting RNA-seq data...")
      melted_rna_marker <<- melt.data.table(
          rnaseq,
          measure.vars = patterns(pattern),
          id.vars = c("gene_id", "gene_name"),
          variable.name = "rna_source",
          value.name = "log2tpm")
      
      print("Melted RNA-seq data.")
    },
    
    merge_rna_meta_and_marker_info = function() {
      melted_rna_marker$gene_id <- gsub("\\..*", "", melted_rna_marker$gene_id)
      if (!is_wide) {
        print("Merging melted RNA-seq data with marker information...")
        
        merged_rna_marker <<- merge.data.table(melted_rna_marker, marker_info,
                                              by = "gene_id", sort = FALSE)
      } else {
        merged_rna_marker <- merge.data.table(melted_rna_marker, marker_info,
                                              by = "gene_id", sort = FALSE)
        merged_rna_marker <<- merge.data.table(merged_rna_marker, rnaseq_meta,
                                              by.x = 'rna_source', by.y = 'sample')
      }
      print(paste("Merged RNA-seq and marker info data with", nrow(merged_rna_marker), "rows."))
    },
      
    calculate_counts = function() {
      print("Calculating counts of genes marked by in different cell lines...")
      counts <<- as.data.table(table(marker_info$cco_levels))
      colnames(counts) <<- c("cco_levels", "COUNTS")
      counts[, cco_levels := as.numeric(cco_levels)]
      print("Counts calculation completed.")
    },
    
    perform_dge_on_first_last_days = function() {
      if (is_wide) {
        print("Skipping DGE analysis for wide data.")
        return()
      }
      
      print("Starting DGE analysis on first and last day samples using counts...")

      # Step 1: Identify the first and last day from the metadata within merged_rna_marker
      rnaseq <<- rnaseq %>% mutate(day_num = as.numeric(gsub("day", "", day)))  # Assuming day column is like 'day1', 'day14', etc.
      
      first_day <- min(rnaseq$day_num)
      last_day <- max(rnaseq$day_num)
      
      print(paste("First day:", first_day))
      print(paste("Last day:", last_day))
      
      # Step 2: Filter for samples from the first and last day
      filtered_data <- rnaseq %>% filter(day_num %in% c(first_day, last_day))
  
      # Step 3: Pivot the data so that the counts for each sample are separate columns
      wide_data <- filtered_data[, c('gene_id', 'count', 'sample')] %>%
        pivot_wider(names_from = sample, values_from = count)
      
      duplicated_columns <- duplicated(wide_data)
      
      if (any(duplicated_columns)) {
        print("Duplicated columns found:")
        print(colnames(wide_data)[duplicated_columns])
        
        # Drop duplicated columns
        wide_data <- wide_data[, !duplicated_columns]
        
        print("Duplicated columns removed.")
      } else {
        print("No duplicated columns found.")
      }

  
      # Step 4: Prepare the DGEList object
      numeric_columns <- setdiff(colnames(wide_data), c("gene_id", "gene_name", "day", "cell_line", "protocol", "day_num", "dgelist_group", "lib_size", "norm_factors"))
      
      dge_counts <- as.matrix(wide_data[, numeric_columns])

      # Create grouping based on day (first_day vs last_day)
      sample_day_map <- filtered_data %>% distinct(sample, day_num) %>% pull(day_num, sample)
      group <- factor(ifelse(sample_day_map[numeric_columns] == first_day, "first_day", "last_day"))
      
      # Step 5: Create the DGEList object
      dge_obj <- DGEList(counts = dge_counts, group = group, genes = wide_data$gene_id)
  
      # Step 6: Normalize and estimate dispersion
      dge_obj <- calcNormFactors(dge_obj)
      design <- model.matrix(~ group)
      dge_obj <- estimateDisp(dge_obj, design)
  
      # Step 7: Fit the model and perform differential expression analysis
      fit <- glmQLFit(dge_obj, design)
      result <- glmQLFTest(fit, coef = 2)
  
      # Step 8: Get the top differentially expressed genes
      top_tags <- topTags(result, n = Inf)
  
      print("DGE analysis complete between first and last day.")
      .self$generate_volcano_plot(top_tags)
    },

    create_ridge_plot = function() {
      message("Generating ridge plot...")
      merged_rna_marker <- na.omit(merged_rna_marker)
      
      log2tpm_col <- ifelse("log2tpm" %in% colnames(merged_rna_marker), "log2tpm", "log2_tpm")
      fill_variable <- ifelse("cell" %in% colnames(merged_rna_marker), "cell", "day")
      merged_rna_marker[[fill_variable]] <- as.factor(merged_rna_marker[[fill_variable]])
      rna_col <- ifelse("rna_source" %in% colnames(merged_rna_marker), "rna_source", "sample")

      unique_rna_sources <- rev(unique(merged_rna_marker[[rna_col]][order(merged_rna_marker[[fill_variable]])]))
      merged_rna_marker[[rna_col]] <- factor(merged_rna_marker[[rna_col]], levels = unique_rna_sources)

      number_of_files_used <- ncol(marker_info[, -c('gene_id')])
      new_col_names <- c()
      for (i in 0:number_of_files_used) {
        new_col_names <- c(new_col_names, setNames(paste('CCO', i), as.character(i)))
      }

      # Start PDF output
      pdf(file = glue("{config$paths$figure_output}/ridge_{config$variables$output_file_name}.pdf"), width = 17, height = 10)

        # Create the ridge plot with the reordered data
        p <- ggplot(data = merged_rna_marker) +
          geom_density_ridges(
            aes_string(x = log2tpm_col, y = rna_col, fill = fill_variable),
            color = "black"
          ) +
          facet_wrap(vars(cco_levels),
            labeller = labeller(cco_levels = as_labeller(new_col_names)),
            nrow = 2
          ) +
          theme(
            axis.text.x = element_text(size = 14, angle = 90),
            axis.title = element_text(size = 15, face = "bold"),
            strip.text.x = element_text(size = 10, face = "bold"),
            plot.title = element_text(size = 22, face = "bold"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15, face = "bold")
          ) +
          scale_fill_manual(name = fill_variable, values = scales::hue_pal()(length(unique(merged_rna_marker[[fill_variable]])))) +
          geom_label(data = counts, aes(x = 12, y = 11, label = COUNTS)) +
          ggtitle(glue("Ridge-plots of gene expression sorted by CCO levels using data from {config$variables$marker_type}"))

        print(p)
        dev.off()
        message("Ridge plot saved to PDF.")
      },

    generate_volcano_plot = function(top_tags) {
      print("Generating volcano plot and merging results...")
      essential_genes <- read.csv(config$paths$essential_genes_file, stringsAsFactors = FALSE)
  
      
      # Merge top_tags with merged_rna_marker
      top_tags_df <- as.data.frame(top_tags)
      top_tags_df$genes <- gsub("\\..*", "", top_tags_df$genes)  # Remove suffix from gene IDs
      
      # Ensure 'gene_id' is consistent in both tables for the merge
      merged_data <- merged_rna_marker[, c('gene_id', 'gene_name', 'cco_levels')]
      merged_data$gene_id <- as.character(merged_data$gene_id)
      
      # Perform the merge operation
      results <- na.omit(merge(
        top_tags_df, merged_data, 
        by.x = 'genes', by.y = 'gene_id', all.x = TRUE
      ))

      # Ensure consistent data types for comparison with essential genes
      results$genes <- as.character(results$genes)
      print(head(results$genes))
      # Debugging: Check if essential genes are properly matched
      results$essential_genes <- results$genes %in% essential_genes$x  # Mark essential genes
      print(table(results$essential_genes))  # Check distribution of TRUE/FALSE values

      # Add other metadata
      results$rank <- rank(results$PValue)
      results$state <- ifelse(results$FDR < 0.05, "sig", "not_sig")
      results$label <- ifelse(results$rank <= 25, results$gene_name, NA_character_)

      # Save results to Excel
      write_xlsx(list(DEGs = results), paste0(config$paths$figure_output, "/deg_results.xlsx"))
      results$essential_genes_label <- ifelse(results$essential_genes, 
                                              "ESSENTIAL GENES", "NON-ESSENTIAL GENES")

      # Generate the Volcano Plot
      pdf(paste0(config$paths$figure_output, "/volcano_plot_essential.pdf"), width = 28, height = 12)

      p <- ggplot(results, aes(x = logFC, y = -log10(PValue), col = -log10(FDR), label = label)) +
        geom_point(size = 4) +
        facet_grid(essential_genes_label ~ cco_levels, scales = "free") +  # Use new labels
        scale_color_gradient2(low = "grey", mid = "grey", high = "red", midpoint = -log10(0.05)) +
        geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
        geom_hline(yintercept = 0, linetype = "dotted", size = 1) +
        geom_text_repel(size = 6, max.overlaps = 15, col = "black") +  # Manage label overlap
        theme_bw() +
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 30, face = "bold"),
          strip.text = element_text(size = 15, face = "bold"),  # Adjust facet label size and style
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 30, face = "bold")
        ) +
        ggtitle("Volcano Plot of DGE with Facets for Essential Genes")

      print(p)
      dev.off()

      print("Volcano plot saved with facets.")

      # Generate line plot for DGE and essentiality
      .self$gen_dge_essentiality_lineplot(results)
    },
    scatter_density_plot = function() {
      plot_with_correlations <- function(plotting_data, first_sample, last_sample) {
        correlation_results <- plotting_data %>%
          group_by(essential_genes, cco_levels) %>%
          summarize(
            spearman = cor(!!sym(first_sample), !!sym(last_sample), method = "spearman", use = "complete.obs"),
            pearson = cor(!!sym(first_sample), !!sym(last_sample), method = "pearson", use = "complete.obs"),
            .groups = 'drop'
          )

        p <- ggplot(plotting_data, aes_string(x = first_sample, y = last_sample)) +
          geom_bin2d(binwidth = c(0.5, 0.5)) +
          scale_fill_viridis_c(limits = c(0, 100)) +
          facet_grid(essential_genes ~ cco_levels, scales = "free") +
          geom_text(data = correlation_results, 
                    aes(x = Inf, y = Inf, label = round(spearman, 2)),
                    inherit.aes = FALSE, vjust = 1, hjust = 1, size = 6, color = "red") +
          theme_minimal() +
          labs(
            x = glue::glue("Log2TPM ({first_sample})"),
            y = glue::glue("Log2TPM ({last_sample})"),
            title = glue::glue("Scatter Plot: {first_sample} vs {last_sample}")
          ) +
          theme(
            strip.text = element_text(size = 15),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 22, face = "bold"),
            plot.title = element_text(size = 30, face = "bold"),
            legend.text = element_text(size = 28),
            legend.title = element_text(size = 30, face = "bold")
          )

        return(list(plot = p, correlation_results = correlation_results))
      }

      essential_genes <- read.csv(config$paths$essential_genes_file, stringsAsFactors = FALSE)

      all_correlation_results <- data.frame()
      merged_rna_marker$day <- as.numeric(gsub("day", "", merged_rna_marker$day))

      first_day <- min(merged_rna_marker$day)
      last_day <- max(merged_rna_marker$day)

      log2tpm_col <- ifelse("log2tpm" %in% colnames(merged_rna_marker), "log2tpm", "log2_tpm")

      first_day_samples <- merged_rna_marker %>%
        filter(day == first_day) %>%
        pull(sample) %>%
        unique()

      last_day_samples <- merged_rna_marker %>%
        filter(day == last_day) %>%
        pull(sample) %>%
        unique()

      plotting_data <- merged_rna_marker %>%
        filter(day %in% c(first_day, last_day)) %>%
        mutate(essential_genes = gene_id %in% essential_genes$x) %>%
        dplyr::select(gene_id, sample, !!log2tpm_col, cco_levels, essential_genes) %>%
        pivot_wider(names_from = sample, values_from = !!log2tpm_col) %>%
        unnest(cols = where(is.list))

      sample_combinations <- expand.grid(
        first_sample = first_day_samples, 
        last_sample = last_day_samples, 
        stringsAsFactors = FALSE
      )

      pdf(paste0(config$paths$figure_output, "/scatter_density_plot_essential.pdf"), width = 24, height = 12)

      for (i in 1:nrow(sample_combinations)) {
        first_sample <- sample_combinations$first_sample[i]
        last_sample <- sample_combinations$last_sample[i]

        if (!(first_sample %in% colnames(plotting_data)) || 
            !(last_sample %in% colnames(plotting_data))) {
          warning(glue::glue("Sample {first_sample} or {last_sample} not found. Skipping."))
          next
        }

        result <- plot_with_correlations(plotting_data, first_sample, last_sample)
        print(result$plot)

        if (!is.null(result$correlation_results)) {
          result$correlation_results <- result$correlation_results %>%
            pivot_longer(cols = c("spearman", "pearson"), 
                        names_to = "correlation_type", 
                        values_to = "cor_value") %>%
            mutate(Comparison = glue::glue("{first_sample} vs {last_sample}"))

          all_correlation_results <- bind_rows(all_correlation_results, result$correlation_results)
        }
      }

      dev.off()

      write.csv(all_correlation_results, 
                paste0(config$paths$figure_output, "/pearson_spearman_correlation_results.csv"), 
                row.names = FALSE)

      print("Scatter density plots saved.")

      .self$generate_spearman_lineplot(all_correlation_results)
    },

    generate_spearman_lineplot = function(correlation_results) {
      print("Generating Spearman and Pearson correlation line plots...")

      # Calculate summary statistics (mean, min, max) for each CCO level and correlation type
      summary_data <- correlation_results %>%
        group_by(correlation_type, essential_genes, cco_levels) %>%
        summarize(
          mean_cor = mean(cor_value, na.rm = TRUE),  # Mean correlation value
          min_cor = min(cor_value, na.rm = TRUE),    # Min correlation value
          max_cor = max(cor_value, na.rm = TRUE),    # Max correlation value
          .groups = 'drop'
        )

      # Create a PDF to store both Pearson and Spearman line plots
      pdf(paste0(config$paths$figure_output, "/corr_lineplot.pdf"), width = 15, height = 15)

      # Loop through each correlation type (Spearman and Pearson) to generate separate plots
      for (cor_type in unique(summary_data$correlation_type)) {
        p <- ggplot(summary_data %>% filter(correlation_type == cor_type),
                  aes(x = cco_levels, y = mean_cor, color = essential_genes)) +
                  geom_line(size = 1.5) +
                  geom_point(size = 3) +
                  geom_ribbon(aes(ymin = min_cor, ymax = max_cor, fill = essential_genes), alpha = 0.2) +
                  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgreen"), 
                                   labels = c("TRUE" = "Essential Genes", "FALSE" = "Non-Essential Genes")) +
                  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkgreen"), 
                                  labels = c("TRUE" = "Essential Genes", "FALSE" = "Non-Essential Genes")) +
                  labs(
                    title = glue("Correlation Distribution ({cor_type}) Across CCO Levels"),
                    x = "CCO Levels", y = glue("{cor_type} Correlation"),
                    color = "Essential Gene", fill = "Essential Gene"
                  ) +
                  theme_minimal() +
                  theme(
                    axis.text = element_text(size = 30),
                    axis.title = element_text(size = 30, face = "bold"),
                    plot.title = element_text(size = 30, face = "bold"),
                    legend.text = element_text(size = 30),
                    legend.title = element_text(size = 30, face = "bold"), 
                    legend.position = c(0.75, 0.1),  # Adjust this to position the legend inside the plot
                    legend.background = element_rect(fill = alpha('white', 0.5)), # Semi-transparent background
                    legend.box.background = element_rect(color = "black", size = 0.5) # Optional border
                  )

        print(p)
      }

    dev.off()
      print("Spearman and Pearson line plots saved.")
    },
    plot_histogram = function() {
      print("Generating histogram of gene counts...")
      pdf(file = glue("{config$paths$figure_output}/{config$variables$output_file_name}.pdf"), width = 17, height = 10)
      print(
        ggplot(counts) + 
          geom_col(aes(x = cco_levels, y = COUNTS), fill = "Orange", color = "black") + 
          xlab(glue('Genes marked by "{config$variables$marker_type}" in X cells')) +
          ylab("Counts") +
          ggtitle(glue('Histogram of counts of genes marked by "{config$variables$marker_type}" in different cell lines')) +
          theme(
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 20, face = "bold"),
            plot.title = element_text(size = 22, face = "bold")
          )
      )
      dev.off()
      print("Histogram saved to PDF.")
    },
    gen_dge_essentiality_lineplot = function(dge_results) {
      # Load DGE Results and Essential Genes
      essential_genes <- read.csv(config$paths$essential_genes_file, stringsAsFactors = FALSE)
  
      
      # Add 'essential_genes' column to DGE results
      dge_results$essential_genes <- dge_results$genes %in% essential_genes$x
      
      # Add FDR significance column
      dge_results$FDR_significant <- dge_results$FDR < 0.05

      # Calculate total gene counts for each CCO level
      total_count <- dge_results %>%
        group_by(cco_levels) %>%
        summarize(total_gene_count = n())
      
      # Prepare lineplot data: Count genes per basket, essentiality, and significance
      lineplot_data <- dge_results %>%
        group_by(cco_levels, essential_genes, FDR_significant) %>%
        summarize(gene_count = n(), .groups = 'drop') %>%
        left_join(total_count, by = "cco_levels") %>%
        mutate(percentage = 100 * (gene_count / total_gene_count))
      # Create PDF to store the plots
      pdf(file = glue("{config$paths$figure_output}/dge_essentials_lineplots.pdf"), width = 15, height = 15)

     p1 <- ggplot(lineplot_data, aes(x = cco_levels, y = gene_count, 
                                color = essential_genes, linetype = FDR_significant)) +
        geom_line(size = 1.5) +  # Thicker lines
        geom_point(size = 4) +   # Larger points
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgreen"),
                          labels = c("TRUE" = "Essential Genes", "FALSE" = "Non-Essential Genes")) +
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted"),
                              labels = c("TRUE" = "Differentially expressed (FDR < 0.05)", "FALSE" = "Not differentially expressed (FDR > 0.05)")) +
        labs(
          title = "Gene Counts Across CCO Levels",
          x = "CCO Levels", y = "Number of Genes",
          color = "Gene Type", linetype = "FDR Significance"
        ) +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 30, face = 'bold'),
          axis.title = element_text(size = 30, face = "bold"),
          plot.title = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30, face = "bold"), 
          legend.position = c(0.5, 0.8),  # Adjust this to position the legend inside the plot
          legend.background = element_rect(fill = alpha('white', 0.5)), # Semi-transparent background
        )

      print(p1)

      # Generate the second plot: Percentage of Genes Across CCO Levels
      p2 <- ggplot(lineplot_data, aes(x = cco_levels, y = percentage, 
                                      color = essential_genes, linetype = FDR_significant)) +
        geom_line(size = 1.5) +  # Thicker lines
        geom_point(size = 4) +   # Larger points
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgreen"),
                          labels = c("TRUE" = "Essential Genes", "FALSE" = "Non-Essential Genes")) +
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted"),
                              labels = c("TRUE" = "Differentially expressed (FDR < 0.05)", "FALSE" = "Not differentially expressed (FDR > 0.05)")) +
        labs(
          title = "Percentage of Genes Across CCO Levels",
          x = "CCO Levels", y = "Percentage of Genes",
          color = "Gene Type", linetype = "FDR Significance"
        ) +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 30, face = 'bold'),
          axis.title = element_text(size = 30, face = "bold"),
          plot.title = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30, face = "bold"),
          legend.position = c(0.5, 0.8),  # Adjust this to position the legend inside the plot
          legend.background = element_rect(fill = alpha('white', 0.5)), # Semi-transparent background
        )

      print(p2)

      # Close the PDF device
      dev.off()

      print("DGE essentiality line plots saved successfully.")
    },
    
    save_results = function() {
      print("Saving melted RNA-seq data and merged RNA-marker data to CSV files...")
      write.csv(melted_rna_marker, glue("{config$paths$figure_output}/melted_table.csv"))
      write.csv(merged_rna_marker, glue("{config$paths$figure_output}/merged_table.csv"))
      print("Results saved successfully.")
    }
  )
)

# Main function to run the analysis
run_analysis <- function(config_file) {
  print("Initializing analysis...")
  analysis <- GeneExpressionAnalysis$new(config_file)
  print("Loading data...")
  analysis$load_data()
  
  print("Calculating cco_levels of true values...")
  analysis$calculate_cco_levels()
  
  print("Melting RNA-seq data...")
  analysis$melt_rna()

  # Always merge and calculate counts
  print("Merging RNA-seq and marker metadata...")
  analysis$merge_rna_meta_and_marker_info()
  
  print("Calculating counts of marked genes...")
  analysis$calculate_counts()
  
  if (!analysis$is_wide) {
    print("Performing DGE analysis on first and last day samples...")
    analysis$perform_dge_on_first_last_days()
    print('Plotting scatter relationship on first and last day...')
    analysis$scatter_density_plot()
  }
  

  print("Creating ridge plots...")
  analysis$create_ridge_plot()
  

  print("Plotting histogram...")
  analysis$plot_histogram()

  print("Saving final results...")
  analysis$save_results()
  
  print("Analysis completed successfully.")
}

# Argument parser
parser <- ArgumentParser()
parser$add_argument("config_file", help = "Path to YAML config file for this script.")
params <- parser$parse_args()

# Run the analysis
run_analysis(params$config_file)