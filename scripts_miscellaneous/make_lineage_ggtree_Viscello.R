##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Oct 24 10:06:42 2020
##########################################################################
##########################################################################
process.bwm.tree.from.viscello = function()
{
  library("treeio")
  library("ggtree")
  
  ##########################################
  # test Viscello code
  ##########################################
  ct_tbl <-  readRDS("../VisCello.celegans/inst/app/data/s6_tbl.rds")
  lin_tbl <-  readRDS("../VisCello.celegans/inst/app/data/s7_tbl.rds")
  #tbl = readRDS("../VisCello.celegans/inst/app/data/lineage_tree_tbl.rds")
  #tree_tbl <- as_tibble(readRDS("../VisCello.celegans/inst/app/data/lineage_tree_tbl.rds"))
  lin_sc_expr <- readRDS("../VisCello.celegans/inst/app/data/lin_sc_expr_190602.rds")
  lin.expanded.list <- readRDS("../VisCello.celegans/inst/app/data/lin_expanded_list_0602.rds")
  avail_nodes <- readRDS("../VisCello.celegans/inst/app/data/avail_nodes.rds")
  
  tree_tbl <- as_tibble(readRDS("../VisCello.celegans/inst/app/data/lineage_tree_tbl.rds"))
  
  ms_tree = tree_tbl;
  root = "MS"; time.cut=400; xmax = 400; color.annot = "lineage"; branch.length='lifetime';
  edge.size = 1; node.lab = TRUE; node.lab.size = 4; node.lab.angle = 90; node.lab.timecut = 80;
  ms_tree<-ms_tree %>% filter(br_time <=time.cut)
  ms_tree <- ms_tree %>% filter(to %in% find_child_root(root, ms_tree))
  ms_tree = as.data.frame(ms_tree)
  
  # remove nodes not to show
  kk = grep('MSa|MSpaa|MSpapa', ms_tree$from)
  jj = grep('MSa|MSpaa|MSpapa', ms_tree$to)
  jj_sel = setdiff(c(1:nrow(ms_tree)), c(kk, jj))
  ms_tree = ms_tree[jj_sel, ]
  
  # change the node names and change the root as well
  ms_tree$from = gsub('MSp', 'MSx', ms_tree$from)
  ms_tree$to = gsub('MSp', 'MSx', ms_tree$to)
  root = 'MS'
  
  kk = which(ms_tree$from == 'MSxapppp')
  ms_tree$to[kk[1]] = 'MSxappppx'
  
  kk = which(ms_tree$from == 'MSxapppa')
  ms_tree$to[kk[1]] = 'MSxapppax'
  
  kk = which(ms_tree$from == 'MSxappa')
  ms_tree$from[kk[1]] = 'MSpappa'
  ms_tree$to[kk[1]] = 'MSpappax'
  
  kk = which(ms_tree$to == 'MSxappa')
  ms_tree$to[kk[1]] = 'MSpappa'
  
  # remove again non-bwm cells
  bwms.all = c('MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
               'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 
               'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
               'MSxpppp', 'MSxpppa', 'MSxppap', #'MSxppaa', 
               'MSxpapp', 'MSxpapa',  
               'MSapaap', 'MSppaap', #'MSxpaap', 
               'MSxpaaa', 'MSxappp', 'MSpappa',  
               'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp', 
               'MSxpp', 'MSxpa', 'MSxap',
               'MSxp', 'MSxa',
               'MSx', 'MS')
  bwm.terminal = c('MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
                   'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 
                   'MSxappppx', 'MSxapppax', 'MSpappax' # all terminal cells
  )
  ms_tree = ms_tree[!is.na(match(ms_tree$from, bwms.all)), ]
  ms_tree = ms_tree[!is.na(match(ms_tree$to, bwms.all)), ]
  
  
  ms_tree$ids = ms_tree$lineage
  ms_tree$ids[!is.na(match(ms_tree$ids, bwm.terminal))] = NA
  
  ms_tree$lineage = ms_tree$to
  ms_tree$value = NA
  
  bwm_tree = as.tibble(ms_tree)
  
  saveRDS(bwm_tree, file = paste0(RdataDir, 'BWM_tree_for_visualization.rds'))
  
}

# functions Viscello
get_tree_color <- function(ctree, pal) {
  tcolor <- colorRampPalette(pal)(length(as_tibble(ctree)$label)) ## (n)
  names(tcolor) <- sort(as_tibble(ctree)$label)
  return(tcolor)
}

make_lineage_ggtree <- function(in_tree = NULL, root = "MS", time.cut = 400, 
                                color.annot = "value", branch.length='lifetime', 
                                tip.lab = F, tip.lab.size = 2, tip.lab.align = F,tip.lab.angle = 0,
                                node.lab = F, node.lab.size = 4, node.lab.timecut = 80, node.lab.angle = 90,
                                xmax = 450, xlim = NULL,
                                color.pal = NULL, edge.size = 1) {
  
  in_tree<-in_tree %>% filter(br_time <=time.cut)
  plot_tree <- in_tree %>% filter(to %in% find_child_root(root, in_tree))
  correct_idx <- which(plot_tree$d_time > xmax)
  plot_tree$d_time[correct_idx] = xmax
  plot_tree$lifetime[correct_idx] <- (plot_tree$d_time - plot_tree$br_time)[correct_idx]
  plot_tree <- as.treedata(plot_tree)
  out_tree <- ggtree(plot_tree, branch.length=branch.length, aes_string(color = color.annot), size = edge.size, ladderize = F) +
    scale_color_continuous(low="black", high="red")
  out_tree <- lineage_tree_flip(out_tree, silence = T)
  
  out_tree = out_tree + geom_tiplab(size = 5) + xlim(0, 380) + 
    geom_label(aes(x=branch, label= ids), fill='gray100', size = 4, nudge_y = 0.3) 

  # if(!is.null(color.pal)) {
  #   out_tree<-out_tree + scale_color_manual(values = color.pal, na.value= "grey")
  # }
  # if(tip.lab) {
  #   out_tree<-out_tree + 
  #     geom_tiplab(size = tip.lab.size,color = "black", align = tip.lab.align, angle = tip.lab.angle)
  # }
  # if(node.lab) {
  #   out_tree<-out_tree + geom_text(aes(x=branch, label=label), data = out_tree$data[out_tree$data$br_time < node.lab.timecut,], 
  #                                  size = node.lab.size, angle=node.lab.angle, color = "black") 
  # }
  # if(!is.null(xlim)) {
  #   out_tree<-out_tree + 
  #     scale_x_continuous(lim = xlim)
  # }
  
  return(out_tree)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


# Computationally check if daugters of all parents have the right p-a order
lineage_tree_flip <- function(treeplot, silence = T) {
  p1 <- treeplot
  for(p in unique(p1$data$parent)) {
    cur_d <- which(p1$data$parent == p)
    cur_d_name <- p1$data$label[cur_d]
    cur_d_n <- p1$data$node[cur_d]
    if(length(cur_d) == 2) {
      cur_d_end <- substrRight(cur_d_name,1)
      cur_d_y <- p1$data$y[cur_d]
      names(cur_d_y) <- cur_d_end
      end_paste <- paste0(cur_d_end, collapse = ",")
      if(end_paste%in% c("a,p", "p,a", "d,v", "v,d", "l,r", "r,l")) {
        if(end_paste%in% c("a,p", "p,a")){
          if(cur_d_y["p"] < cur_d_y["a"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(end_paste%in% c("d,v", "v,d")){
          if(cur_d_y["v"] < cur_d_y["d"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(end_paste%in% c("l,r", "r,l")){
          if(cur_d_y["r"] < cur_d_y["l"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
      } else {
        cur_d_y <- p1$data$y[cur_d]
        names(cur_d_y) <- cur_d_name
        if(all(c("Z2", "Z3") %in% cur_d_name)) {
          if(cur_d_y["Z2"] < cur_d_y["Z3"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(all(c("P2", "EMS") %in% cur_d_name)) {
          if(cur_d_y["P2"] < cur_d_y["EMS"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(all(c("E", "MS") %in% cur_d_name)) {
          if(cur_d_y["E"] < cur_d_y["MS"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(all(c("C","P3") %in% cur_d_name)) {
          if(cur_d_y["P3"] < cur_d_y["C"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
        if(all(c("D", "P4") %in% cur_d_name)) {
          if(cur_d_y["P4"] < cur_d_y["D"]) {
            if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
            p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
          }
        }
      }
    }
  }
  return(p1)
}


find_child_root <- function(root, tree_tbl) {
  find_child <- function(root, tree_tbl) {
    dts <- tree_tbl$to[which(tree_tbl$from == root)]
    return(c(dts, lapply(dts, find_child, tree_tbl)))
  }
  
  res<-find_child(root, tree_tbl)
  res <- unlist(res, recursive = T)
  return(res)
}


# Function to get lowest common ancester
get_ancestor <- function(df, node){
  anc <- get_parent(df, node)
  i <- 1
  while (i <= length(anc)) {
    anc <- c(anc, get_parent(df, anc[i]))
    i <- i + 1
  }
  return(anc)
}


get_parent <- function (df, node) 
{
  parent_id <- df$from[df$to == node]
  parent_id[parent_id != node]
}


get_level <- function (df, node, root = "P0") 
{
  lev = 0
  cur_n <- node
  parent  <- node
  while(parent != root) {
    parent <- get_parent(df, cur_n)
    cur_n <- parent
    lev <- lev+1
  }
  return(lev)
}


