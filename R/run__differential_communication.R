
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
input_path_ligand_receptor_count_output = args[2:(length(args)-5)]
input_metadata = args[[length(args)-4]]
output_path_plot_overall = args[[length(args)-3]]
output_path_plot_heatmap = args[[length(args)-2]]
output_path_plot_circle = args[[length(args)-1]]
output_path_values_communication = args[[length(args)]]

renv::load(project = code_directory)

library(tidyverse)
library(Seurat)
library(glue)
library(CellChat)
library(tidyseurat)
library(tidySingleCellExperiment)
library(tidyHeatmap)
library(purrr)
library(patchwork)
library(cowplot)
library(gridGraphics)

source(glue("{code_directory}/R/in_house_functions.R"))

# Create dir
output_path_plot_overall |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

factor_of_interest = as.symbol("factor_of_interest")

# input_path_ligand_receptor_count_output = dir("data/fast_pipeline_results/communication", full.names = TRUE)

integrated_cellchat_plot_df =
  input_path_ligand_receptor_count_output |>
  map_df(~ readRDS(.x)) |>

  # Add metadata
  left_join( readRDS("data/metadata.rds"), by = "sample") |>

  # Nest
  nest(data_DB = -DB) |>
  mutate(joint = map(data_DB, ~ {
    object.list <- .x$data
    mergeCellChat(object.list, add.names = paste(.x$sample, pull(.x, !!factor_of_interest), sep="___"))
  })) |>

  # Add histogram
  mutate(plot_histogram = map2(joint, DB, ~ {
   s =  sapply(.x@net, function(x) sum(x$count))

    when(is.na(s) ~ 0 , ~s) |>
    as.data.frame() |>
    as_tibble(rownames = quo_name(factor_of_interest)) |>
    setNames(c( quo_name(factor_of_interest), "tot"))
  })) |>
  select(-joint) |>

  # Information flow
  # Add histogram
  mutate(plot_information_flow = map2(plot_histogram,DB, ~ {

    .x	|>
      separate(!!factor_of_interest, c(  "sample", quo_name(factor_of_interest)), "___") |>
      ggplot(aes(!!factor_of_interest, tot, fill=!!factor_of_interest)) +
      geom_boxplot() +
      geom_point() +
      ggpubr::stat_compare_means() +
      ggtitle(.y)
  }
  ))

# Plot overall counts
p =
  integrated_cellchat_plot_df |>
  pull(plot_information_flow) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = 'collect' ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

ggsave(
  plot = p,
  filename = output_path_plot_overall,
  units = "mm",
  width = 60*sqrt(length(p)),
  height = 60*sqrt(length(p)),
  device = "pdf"
)


# integrated_cellchat_plot_df |> saveRDS("cancer_only_analyses/integrated_cellchat_plot_df_sample_wise.rds")

values_df_for_heatmap =
  integrated_cellchat_plot_df |>
  select(DB, data_DB) |>
  unnest(data_DB) |>
  mutate(data = map2(
    data  , sample,
    ~ when(
      .x,
      length(.x@netP$pathways) > 0 ~
        netAnalysis_computeCentrality(., slot.name = "netP") |>
        cellchat_process_sample_signal(
          pattern = "all", signaling = .x@netP$pathways,
          title = .y, width = 5, height = 6, color.heatmap = "OrRd"
        ),
      ~ (.)
      )
  )) |>
  unnest(data)  |>
  complete(nesting(DB, gene), nesting(!!factor_of_interest, sample), cell_type, fill = list(value = 0))

# Save for machine learning
values_df_for_heatmap |> saveRDS(output_path_values_communication)

signal_df_for_heatmap =
  values_df_for_heatmap |>
  nest(value_data = -c(DB  ,   gene ,  !!factor_of_interest,  cell_type)) |>
  mutate(value = map_dbl(value_data, ~ mean(.x$value))) |>

  select(-value_data ) |>
  spread(!!factor_of_interest, value) |>

  # THIS SHOULD CHANGE WITH THE RIGHT VALUES OF FACTOR OF INTEREST
  mutate(diff = high - low) |>

  # Fix CD40
  mutate(gene = if_else(gene=="CD40" & DB=="Cell-Cell Contact" , "CD40_",  gene)) |>
  group_by(gene) |>
  mutate(max_diff = max(abs(diff))) |>
  filter(abs(max_diff)>0.34)

# Plot
signal_df_for_heatmap |>
  group_by(cell_type) |>
  mutate(sum_cell_diff = sum(diff)) |>
  group_by(gene) |>
  mutate(sum_gene_diff = sum(diff)) |>
  group_by(DB) |>

   mutate(gene = if_else(gene=="CD40" & DB=="Cell-Cell Contact" , "CD40_",  gene)) |>

  heatmap(
    gene, cell_type, diff,
    palette_value = circlize::colorRamp2(seq(0.5, -0.5, length.out =11), RColorBrewer::brewer.pal(11, "RdBu")),
    scale="none"
  ) |>
  add_bar(sum_cell_diff) |>
  add_bar(sum_gene_diff) |>
  save_pdf(output_path_plot_heatmap)



# plot_circla_communication_all
plot_circla_communication_all =
  integrated_cellchat_plot_df |>

  unnest(data_DB) |>
  left_join(signal_df_for_heatmap |> distinct(DB, gene) |> ungroup()) |>
  mutate(genes_in_pathway = map2(gene, data, ~ select_genes_for_circle_plot(.y, .x))) |>

  mutate(plot_cell_Cell_pathway = map2( data, gene , ~
                                          cellchat_matrix_for_circle(.x,  layout = "circle", signaling = .y)
  ))  |>
  filter(map_lgl(plot_cell_Cell_pathway, ~ !is.null(.x))) |>

  # Create the full column name because cellchat sucks
  mutate(all_cell_types = map(plot_cell_Cell_pathway, ~ colnames(.x))) |>
  mutate(all_cell_types = unlist(all_cell_types) |> unique() |> list()) |>
  mutate(missing_column = map2(
    plot_cell_Cell_pathway,
    all_cell_types,
    ~ .y |> setdiff(colnames(.x))
  )) |>
  mutate(plot_cell_Cell_pathway =
           map2(
             plot_cell_Cell_pathway,
             missing_column,
             ~ when(.y,
                    length(.) == 0 ~ .x,
                    ~ {
                      x = as.data.frame(.x)
                      x[.y] = 0
                      x[.y, ] = 0
                      x = as.matrix(x)

                    })
           )) |>
  mutate(plot_cell_Cell_pathway = map2(plot_cell_Cell_pathway, all_cell_types, ~ .x[.y, .y])) |>


  # Add line weights
  mutate(
    line_weights = map2(
      data,
      all_cell_types,
      ~ .y |>
        purrr::map_dfc(setNames, object = list(integer())) |>
        bind_rows(
          .x@idents |> table() |>  enframe() |> spread(name, value) |> mutate(across(everything(), ~ as.integer(.x)))
        ) |>
        mutate(across(everything(), ~ replace_na(.x, 0)))
    )
  ) |>

  # Summarise
  nest(data_path_type = -c(gene, !!factor_of_interest)) |>
  mutate(
    pathway_mean = map(data_path_type, ~ purrr::reduce(.x$plot_cell_Cell_pathway, `+`)),
    line_weights_sum = map(
      data_path_type,
      ~ .x$line_weights |> purrr::reduce(bind_rows) |> colSums()
    ),
    genes_in_pathway = map(
      data_path_type,
      ~ .x$genes_in_pathway |> unlist() |> unique()
    )
  ) |>
  arrange(!!factor_of_interest) |>
  nest(data_type_path = -c(gene)) |>
  mutate(
    plot_diff = map(
      data_type_path,
      ~  when(
        pull(.x, !!factor_of_interest),
        unlist(.) |>  identical("low") ~ -.x$pathway_mean[[1]],
        unlist(.) |>  identical("high") ~ .x$pathway_mean[[1]],
        unlist(.) |> identical(c("high", "low")) ~ .x$pathway_mean[[1]] - .x$pathway_mean[[2]],
        ~ stop("something went wrong")
      )
    ),
    line_weights_sum_sum = map(
      data_type_path,
      ~ .x$line_weights_sum |> purrr::reduce(`+`) |> magrittr::divide_by(2)
    ),
    genes_in_pathway = map(
      data_type_path,
      ~ .x$genes_in_pathway |> unlist() |> unique()
    )
  ) |>
  mutate(plot_diff_max = map_dbl(plot_diff, ~ max(abs(.x)))) |>
  mutate(plot_diff_quant = quantile(plot_diff_max, 0.25)) |>
  mutate(circle_plot = pmap(
    list(plot_diff, line_weights_sum_sum, gene, genes_in_pathway, plot_diff_quant),
    ~ {print(".");

        draw_cellchat_circle_plot(
        ..1,
        vertex.weight = ..2,
        title.name =  paste(..3, "\n", ..4),
        edge.width.max = 4,
        remove.isolate = TRUE,
        top_absolute=..5,
        top = 0.2,
        arrow.width = 4
      )
    }
  ))

p = plot_circla_communication_all |>
  pull(circle_plot) |>
  discard(is.null)

ggsave(
  plot = purrr::reduce(p, `+`) + plot_layout( ncol = ceiling(sqrt(length(p))), nrow = ceiling(sqrt(length(p)))),
  filename = output_path_plot_circle,
  units = "mm",
  width = 60*sqrt(length(p)),
  height = 60*sqrt(length(p)),
  device = "pdf"
)

