library(tidyverse)
library(here)
library(data.table)

rm(list = ls())

##################################################################################
analysisName <- "cleistothecia_size"
outDir <- here::here("analysis", "00_temp")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_data <- here("analysis", "00_temp", "inputfile_cleistothecia_size.txt")

inputData <- suppressMessages(readr::read_tsv(file = file_data))

plotDf <- data.table::melt(data = as.data.table(inputData), id.vars = "time",
                 variable.name = "gene", value.name = "size") %>% 
  dplyr::as_tibble() %>% 
  dplyr::group_by(gene, time) %>% 
  dplyr::summarise(
    avg_size = mean(size),
    sd = sd(size)
  ) %>% 
  dplyr::ungroup()


barWd <- 0.5

pt <- ggplot(data = plotDf, mapping = aes(x = gene, y = avg_size, fill = time)) +
  geom_bar(stat = "identity", position = position_dodge(width = barWd),
           width = 0.5) +
  geom_errorbar(mapping = aes(ymin = avg_size - sd, ymax = avg_size + sd),
                position = position_dodge(width = barWd)) +
  geom_rect(data = dplyr::filter(plotDf, gene == "WT"),
            mapping = aes(xmin = as.numeric(gene) - 0.5, xmax = as.numeric(gene) + 0.5,
                          ymin = -Inf, ymax = Inf),
            fill = "#a6a6a6", alpha = .2) +
  labs(title = "Cleistothecia size",
       y = "Cleistothecia size (mm)") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
  scale_fill_manual(
    name = "Time",
    values = c("10d" = "#a3cbb8", "14d" = "#4f5d75")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))


png(filename = paste(outPrefix, ".png", sep = ""), width = 4000, height = 2000, res = 400)
pt
dev.off()

pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 14, height = 7)
pt
dev.off()

