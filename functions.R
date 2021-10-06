# functions for nil browser shiny
library(tidyverse)

# plot NILs
nil_plot <- function(data, chr, left.cb = 0, left.n2 = 0, left.bound = 0, right.bound = 19e6, scan.range = 2e4, 
                     all.chr=F, section = "all", background = F, ci = NA, elements = F, order = T, n2_color = "orange", cb_color = "blue",
                     tsize = 12){
  
  # copy nilgeno
  nilgeno2 <- data %>%
    dplyr::mutate(sample = strain) %>%
    dplyr::mutate(end = as.numeric(end),
                  chrom = chromosome,
                  gt_name = nil_genotype)
  strains <- unique(data$strain)
  
  # check to see if all strains are represented in the nilgeno2 file
  df <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::distinct(sample)
  
  if(nrow(df) != length(strains)) {
    unknown_strains <- setdiff(strains, df$sample)
    
    # make fake genotype for unknown strains
    for(s in unknown_strains) {
      fake_geno <- nilgeno2 %>%
        dplyr::filter(sample == "N2")%>%
        dplyr::distinct(chrom, .keep_all = T) %>%
        dplyr::mutate(sample = s,
                      gt = 3,
                      gt_name = "unknown")
      
      # add fake geno to nil geno
      nilgeno2 <- rbind(nilgeno2, fake_geno)
    }
  }
  
  # # # determine if NILs are CB or N2
  nilsII_sort_type <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::group_by(sample)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(sample, gt)%>%
    dplyr::mutate(major_gt = sum(gt_ct))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(major_gt))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::mutate(nil_type = ifelse(gt == 1, "CB", "N2"))%>%
    dplyr::select(sample, nil_type)
  
  # # # keep NILs that lost NIL genotype on right side
  nilsII_left <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::filter((start > left.cb - .3e6 & start < left.cb + .3e6 & nil_type == "CB") |
                    (start > left.n2 - .3e6 & start < left.n2 + .3e6 & nil_type == "N2") )%>%
    dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
    dplyr::mutate(side = "LEFT")%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(size))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::filter(size > scan.range)%>% # # # remove small (likely wrong calls) around interval site
    dplyr::select(sample, side)
  
  # # # keep NILs that lost NIL genotype on left side
  nilsII_right <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::filter(!sample %in% nilsII_left$sample)%>%
    dplyr::mutate(side = "RIGHT")%>%
    dplyr::distinct(sample,.keep_all = T)%>%
    dplyr::select(sample, side)
  
  nil_sides <- bind_rows(nilsII_left,nilsII_right)
  
  
  nilsII_sort_left <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::left_join(.,nil_sides, by = "sample")%>%
    dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
    dplyr::group_by(sample)%>%
    dplyr::filter(start > left.bound & start < right.bound)%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::filter(side == "LEFT")%>%
    dplyr::arrange( desc(gt_ct))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::arrange(nil_type, desc(gt_ct))
  
  #Here
  nilsII_sort_right <- nilgeno2 %>%
    dplyr::filter(sample %in% strains) %>%
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::left_join(.,nil_sides, by = "sample")%>%
    dplyr::filter(side == "RIGHT")%>%
    dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
    dplyr::group_by(sample)%>%
    dplyr::filter(start > left.bound & start < right.bound)%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::arrange( desc(gt_ct))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::ungroup()%>%
    dplyr::arrange( nil_type, gt_ct)
  
  N2CB <- nilgeno2 %>%
    dplyr::filter(sample %in% c("N2", "CB4856"), chrom == chr) %>%
    dplyr::mutate(nil_type = "parent", side = NA, size = NA, gt_ct = NA)
  
  nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left) %>%
    arrange(desc(nil_type), desc(side), size)
  
  nilsII_sort <- rbind(nilsII_sort, N2CB)
  
  if(nrow(df) != length(strains)) {
    nilsII_sort <- rbind(nilsII_sort, nilgeno2 %>%
                           dplyr::filter(sample %in% unknown_strains,
                                         chrom == chr) %>%
                           dplyr::mutate(nil_type = "unkown", side = NA, size = NA, gt_ct = NA))
    
  }
  
  if (all.chr == T){
    
    nilsII <- nilgeno2 %>%
      dplyr::filter(sample %in% strains) %>%
      dplyr::filter(chrom != "MtDNA")%>%
      dplyr::left_join(.,nilsII_sort_type, by = "sample")
    
    if(nrow(df) != length(strains)) {
      nilsII <- nilsII %>%
        dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
    }
    
    if(order == T) {
      nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
    } else {
      nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
    }
    nilsII$gt <- as.character(nilsII$gt)
  } else {
    nilsII <- nilgeno2 %>%
      dplyr::filter(sample %in% strains) %>%
      dplyr::filter(chrom == chr, chrom != "MtDNA")%>%
      dplyr::left_join(.,nilsII_sort_type, by = "sample")
    
    if(nrow(df) != length(strains)) {
      nilsII <- nilsII %>%
        dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
    }
    
    if(order == T) {
      nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
    } else {
      nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
    }
    nilsII$gt <- as.character(nilsII$gt)
  }
  
  # make fake 'background' genome
  if(background == T) {
    # cannot show all chromosomes and the background
    if(all.chr == T) {
      bgplot <- NULL
    } else {
      bg <- data.frame(chrom = NA, start = NA, end = NA, sample = NA, gt = NA, gt_name = NA, supporting_sites = NA, sites = NA,
                       DP = NA, switches = NA, nil_type = NA)
      for(i in unique(nilsII$sample)) {
        # get the NIL type
        type <- nilsII %>%
          dplyr::filter(sample == i) %>%
          dplyr::distinct(nil_type) %>%
          dplyr::pull(nil_type)
        
        # add a new row with the background genotype in a "new" chromosome
        bg <- rbind(bg, c("Genome", 1, 3e6, i, 
                          ifelse(type == "N2", 2,
                                 ifelse(type == "CB", 1,
                                        ifelse(type == "unknown", 3, NA))),
                          ifelse(type == "N2", "CB4856",
                                 ifelse(type == "CB", "N2",
                                        ifelse(type == "unknown", "unknown", NA))),
                          NA, NA, NA, NA, type))
      }
      bg <- bg %>%
        tidyr::drop_na(sample) %>%
        dplyr::mutate(start = as.numeric(start), end = as.numeric(end))
      
      # add "chrom" to chromosome for plotting
      nils2 <- nilsII %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
      
      # background plot
      bgplot <- ggplot(bg)+
        geom_segment(aes(x = start/1e6, y = factor(sample, levels = levels(nils2$sample)), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
        scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
        # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
        facet_grid(~chrom, scales = "free",  space = "free")+
        theme_bw(tsize) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              strip.text = element_text(face = "bold", color = "black"),
              plot.title = element_text(face="bold"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.ticks = element_blank())
    }
  } else { bgplot <- NULL }
  
  # only plot the N2-NILs
  if(section == "N2-NILs") {
    nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
      geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
      facet_grid(~chrom, scales = "free",  space = "free")+
      scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
      theme_bw(tsize) +
      theme(axis.text.x = element_text(face="bold", color="black"),
            axis.text.y = element_text(face="bold", color="black"),
            axis.title.x = element_text(face="bold", color="black"),
            axis.title.y = element_text(face="bold", color="black"),
            strip.text = element_text(face = "bold", color = "black"),
            plot.title = element_text(face="bold"),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      geom_vline(data = filter(nilsII, chrom == chr) %>%
                   dplyr::distinct(sample, chrom) %>%
                   dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                   tidyr::separate_rows(xint), 
                 aes(xintercept = as.numeric(xint)), color = "red", size = 2) +
      labs(x = "Genomic position (Mb)", y = "")
  } else if(section == "CB-NILs") {
    # only plot the CB-NILs
    nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
      geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
      facet_grid(~chrom, scales = "free",  space = "free")+
      scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
      theme_bw(tsize) +
      theme(axis.text.x = element_text(face="bold", color="black"),
            axis.text.y = element_text(face="bold", color="black"),
            axis.title.x = element_text(face="bold", color="black"),
            axis.title.y = element_text(face="bold", color="black"),
            strip.text = element_text(face = "bold", color = "black"),
            plot.title = element_text(face="bold"),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      geom_vline(data = filter(nilsII, chrom == chr) %>%
                   dplyr::distinct(sample, chrom) %>%
                   dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                   tidyr::separate_rows(xint), 
                 aes(xintercept = as.numeric(xint)), color = "red") +
      labs(x = "Genomic position (Mb)", y = "")
  } else {
    
    # plot all NILs
    nl.pl <- ggplot(nilsII)+
      geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
      facet_grid(~chrom, scales = "free",  space = "free")+
      scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
      theme_bw(tsize) +
      theme(axis.text.x = element_text(face="bold", color="black"),
            axis.text.y = element_text(face="bold", color="black"),
            axis.title.x = element_text(face="bold", color="black"),
            axis.title.y = element_text(face="bold", color="black"),
            strip.text = element_text(face = "bold", color = "black"),
            plot.title = element_text(face="bold"),
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      geom_vline(data = filter(nilsII, chrom == chr) %>%
                   dplyr::distinct(sample, chrom) %>%
                   dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                   tidyr::separate_rows(xint), 
                 aes(xintercept = as.numeric(xint)), color = "red", size = 2) +
      labs(x = "Genomic position (Mb)", y = "")
  }
  
  # return plots and dataframe
  if(is.null(bgplot)) {
    return(list(nl.pl, nilsII_sort, nilsII))
  } else {
    if(elements == T) {
      return(list(nl.pl, bgplot, nilsII_sort))
    } else {
      return(list(cowplot::plot_grid(nl.pl, bgplot, nrow = 1, ncol = 2, align = "h", axis = "b", rel_widths = c(1, 0.3)),
                  nilsII_sort,
                  nilsII,
                  nilgeno2))
      # return(list(nl.pl,
      #             nilsII_sort,
      #             nilsII,
      #             bgplot))
    }
  }
}
