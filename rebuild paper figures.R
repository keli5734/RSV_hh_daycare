
library(tidyverse)
library(ggsci)

## ------------------------------ CONFIG ------------------------------------
in_dir  <- "Results"
out_dir <- in_dir
B_BOOT       <- 1000L
PREV_METRIC  <- "sym"                 # "all" incl. asymptomatic I_ab ; "sym" symptomatic
CONTRAST_COVS <- c(0, 100)
SETTINGS_KEEP <- c("daycare","household")
EXCLUDE_STAFF_IN_HH   <- TRUE
MERGE_STAFF_INTO_ADULT   <- FALSE
MERGE_ELDERLY_INTO_ADULT <- TRUE
CONTRAST_MIN_EVENTS   <- 0.10
ROLES <- c("infant","toddler","adult","elderly","staff")
month_starts <- c(1,32,63,93,124,154,184,215,243,274,304,335)
month_labels <- c("Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun")
icol <- c(infant="#e6550d", toddler="#31a354", adult="#3182bd", elderly="#756bb1", staff="#fa9fb5")
have_gg <- requireNamespace("ggplot2",   quietly = TRUE); if (have_gg) library(ggplot2)
have_al <- requireNamespace("ggalluvial", quietly = TRUE)
theme_pub <- function() theme_classic(base_size = 13) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"),
        plot.caption = element_text(size = 7.5, colour = "grey45", hjust = 0)) +
  theme(
    text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 15)
  )
cov_of <- function(f) as.integer(sub(".*cov(\\d{3})\\.rds$", "\\1", basename(f)))
setting_of <- function(s) ifelse(grepl("^daycare", s), "daycare",
                                 ifelse(grepl("^household", s), "household", "community"))
files <- list.files(in_dir, "^RSV_rawsims_cov\\d{3}\\.rds$", full.names = TRUE)
if (length(files) == 0) stop(sprintf("No packs in '%s/'.", in_dir))
files <- files[order(cov_of(files))]; covs <- cov_of(files)
col_inf <- if (PREV_METRIC == "all") "inf_dc_I" else "inf_dc_Isym"
col_tod <- if (PREV_METRIC == "all") "tod_dc_I" else "tod_dc_Isym"
col_adu <- if (PREV_METRIC == "all") "hh_I"     else "hh_I_sym"
y_lab   <- if (PREV_METRIC == "all") "Number of infections" else  "Number of infections"



cnt_list <- list(); curve_sys <- list()
flow_raw <- list(); npacks_by_cov <- integer(0)
for (f in files) {
  cov <- cov_of(f); ck <- as.character(cov)
  packs <- readRDS(f); np <- length(packs)
  nd <- nrow(packs[[1]]$daily)
  sims <- vapply(packs, function(p) if (!is.null(p$sim)) as.integer(p$sim) else NA_integer_, integer(1))


  cn <- matrix(0, np, 21)
  colnames(cn) <- c("cd_vaxinf","n_vaxinf","cd_unvaxinf","n_unvaxinf","cd_tod","n_tod",
                    "cd_adult","n_adult","cd_eld","n_eld",
                    "cd_adulteld","n_adulteld",
                    "ci_vaxinf","ci_unvaxinf","ci_tod","ci_adult","ci_eld","ci_adulteld",
                    "cd_total","n_total","ci_total")
  for (r in seq_len(np)) {
    af <- packs[[r]]$agents_final
    ever <- if (!is.null(af$ever_infected)) af$ever_infected else af$n_infections > 0
    staff_hh <- unique(af$hh_id[af$role == "staff"]); core <- !(af$hh_id %in% staff_hh)
    mvi <- af$intended_vax & af$role=="infant" & af$is_daycare
    mui <- !af$intended_vax & af$role=="infant" & af$is_daycare
    mt  <- af$role=="toddler"; ma <- af$role=="adult" & core; me <- af$role=="elderly" & core
    mae <- ma | me                                              # <-- combined adult+elderly mask
    cn[r,] <- c(sum(af$infected_outcome[mvi]), sum(mvi),
                sum(af$infected_outcome[mui]), sum(mui),
                sum(af$infected_outcome[mt]),  sum(mt),
                sum(af$infected_outcome[ma]),  sum(ma),
                sum(af$infected_outcome[me]),  sum(me),
                sum(af$infected_outcome[mae]), sum(mae),
                sum(ever[mvi]), sum(ever[mui]), sum(ever[mt]), sum(ever[ma]), sum(ever[me]), sum(ever[mae]),
                sum(af$infected_outcome), nrow(af), sum(ever))
  }
  cnt_list[[ck]] <- as.data.frame(cn)

   usims <- sort(unique(sims[!is.na(sims)]))
  grp <- if (length(usims)==0) as.list(seq_len(np)) else lapply(usims, function(s) which(sims==s))
  gsz <- vapply(grp, length, integer(1))
  if (length(unique(gsz)) > 1L)
    warning(sprintf("cov%03d: uneven #centres per sim group (%d-%d).", cov, min(gsz), max(gsz)))
  INF <- TOD <- ADU <- matrix(0, length(grp), nd)
  for (si in seq_along(grp)) {
    ai <- at <- aa <- numeric(nd)
    for (r in grp[[si]]) { d <- packs[[r]]$daily
    ai <- ai + d[[col_inf]]; at <- at + d[[col_tod]]; aa <- aa + d[[col_adu]] }
    INF[si,] <- ai; TOD[si,] <- at; ADU[si,] <- aa
  }
  curve_sys[[ck]] <- list(inf=INF, tod=TOD, adu=ADU)

   if (cov %in% CONTRAST_COVS) {
    npacks_by_cov[ck] <- np
    acc <- list()
    for (r in seq_len(np)) {
      af <- packs[[r]]$agents_final; lg <- packs[[r]]$infection_log
      if (is.null(lg) || nrow(lg) == 0) next
      rmap <- af$role; names(rmap) <- as.character(af$global_id)
      ir <- ifelse(is.na(lg$infector), "community", rmap[as.character(lg$infector)])
      ee <- lg$role_at_infection; st <- setting_of(lg$source)
      if (MERGE_STAFF_INTO_ADULT)   { ir[ir=="staff"]   <- "adult"; ee[ee=="staff"]   <- "adult" }
      if (MERGE_ELDERLY_INTO_ADULT) { ir[ir=="elderly"] <- "adult"; ee[ee=="elderly"] <- "adult" }
      ok <- ee %in% ROLES & ir %in% ROLES & st %in% SETTINGS_KEEP
      if (EXCLUDE_STAFF_IN_HH)
        ok <- ok & !(st == "household" & (ir == "staff" | ee == "staff"))
      if (!any(ok)) next
      key <- paste(ir[ok], ee[ok], st[ok], sep="|"); tb <- table(key)
      for (nm in names(tb)) acc[[nm]] <- (if (is.null(acc[[nm]])) 0L else acc[[nm]]) + tb[[nm]]
    }
    if (length(acc)) {
      p <- do.call(rbind, strsplit(names(acc), "|", fixed=TRUE))
      flow_raw[[ck]] <- data.frame(coverage=cov, infector=p[,1], infectee=p[,2],
                                   setting=p[,3], count=as.integer(unlist(acc)),
                                   stringsAsFactors=FALSE)
    }
  }
  rm(packs); gc(verbose = FALSE)
  cat(sprintf("cov%03d: %d packs, %d sim-groups (centres/group %d-%d)\n",
              cov, np, length(grp), min(gsz), max(gsz)))
}

## ================================ FIGURE 1 =================================
cov_cols <- c(
  "0%"   = "#000000",  # black
  "20%"  = "#0072B2",  # blue
  "40%"  = "#56B4E9",  # sky blue
  "60%"  = "#009E73",  # green
  "80%"  = "#E69F00",  # orange
  "100%" = "#D55E00"   # red-orange
)
if (have_gg) {
  panel_titles <- c(inf="A) Infants", tod="B) Toddlers", adu="C) Adults & older adults")
  nd <- ncol(curve_sys[[1]]$inf)
  f1 <- do.call(rbind, lapply(names(curve_sys), function(k) {
    cs <- curve_sys[[k]]
    do.call(rbind, lapply(c("inf","tod","adu"), function(s) {
      m <- cs[[s]]
      data.frame(day=seq_len(nd), coverage=as.integer(k), panel=s,
                 mean=colMeans(m), lo=apply(m,2,quantile,0.025), hi=apply(m,2,quantile,0.975),
                 row.names=NULL) })) }))
  f1$panel_lab <- factor(panel_titles[f1$panel], levels=panel_titles)
  cvs <- sort(unique(f1$coverage)); f1$cov_lab <- factor(paste0(f1$coverage,"%"), levels=paste0(cvs,"%"))
  g1 <- ggplot(f1, aes(day, mean, colour=cov_lab, fill=cov_lab, group=cov_lab)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.12, colour=NA) +
    geom_line(linewidth=1) +
    facet_wrap(~ panel_lab, ncol=1, scales="free_y",
               labeller = as_labeller(c(
                 "Infants"                 = "A) Infants",
                 "Toddlers"                = "C) Toddlers",
                 "Adults & older adults"   = "E) Adults & older adults"))) +
    scale_x_continuous(breaks=month_starts, labels=month_labels, expand=c(0.01,0)) +
    scale_colour_manual(values = cov_cols, name = "mAb coverage") +
    scale_fill_manual(values = cov_cols, name = "mAb coverage") +
    labs(x="Month", y=y_lab) +
    theme_pub() + theme(strip.text=element_text(hjust=0), legend.position = "right")
  ggsave(file.path(out_dir,"FIG1_prevalence_system.pdf"), g1, width=8.5, height=7)
  write.csv(f1[,c("coverage","panel","day","mean","lo","hi")],
            file.path(out_dir,"FIG1_prevalence_system.csv"), row.names=FALSE)
  cat("[fig1] written\n")
}

## ================================ FIGURE 2 =================================
ARd <- function(cn,g){ t<-sum(cn[[paste0("n_",g)]]); if(t>0) sum(cn[[paste0("cd_",g)]])/t else NA_real_ }
ARi <- function(cn,g){ t<-sum(cn[[paste0("n_",g)]]); if(t>0) sum(cn[[paste0("ci_",g)]])/t else NA_real_ }
base_cn <- cnt_list[[as.character(min(covs))]]
prot_point <- function(cn, base) c(
  direct=1-ARd(cn,"vaxinf")/ARd(cn,"unvaxinf"), total=1-ARd(cn,"vaxinf")/ARd(base,"unvaxinf"),
  tod=1-ARd(cn,"tod")/ARd(base,"tod"),
  adulteld=1-ARd(cn,"adulteld")/ARd(base,"adulteld"))                 # <-- combined, was adult= + eld=
boot_prot <- function(cn, base, B=B_BOOT) {
  est <- prot_point(cn, base); M <- matrix(NA_real_, B, length(est))
  for (b in seq_len(B)) {
    cb <- cn[sample.int(nrow(cn), replace=TRUE),,drop=FALSE]
    bb <- base[sample.int(nrow(base), replace=TRUE),,drop=FALSE]
    M[b,] <- prot_point(cb, bb) }
  data.frame(metric=names(est), est=est, lo=apply(M,2,quantile,0.025,na.rm=TRUE),
             hi=apply(M,2,quantile,0.975,na.rm=TRUE), row.names=NULL)}
prot_df <- do.call(rbind, lapply(covs, function(cv){
  d <- boot_prot(cnt_list[[as.character(cv)]], base_cn); d$coverage <- cv; d }))
lab <- c(direct="Direct (infant)", total="Total (infant)", tod="Toddler",
         adulteld="Adult & older adult")                              # <-- one combined label
pnl <- c(direct="A) Infant (disease)", total="A) Infant (disease)",
         tod="B) Indirect (disease)", adulteld="B) Indirect (disease)")
prot_df$label <- lab[prot_df$metric]; prot_df$panel <- pnl[prot_df$metric]
write.csv(prot_df, file.path(out_dir,"FIG2_protection.csv"), row.names=FALSE)
if (have_gg) {
  pp <- prot_df[!is.na(prot_df$est), ]; pp[,c("est","lo","hi")] <- pp[,c("est","lo","hi")]*100
  pp$panel <- factor(pp$panel, levels=c("A) Infant (disease)","B) Indirect (disease)"))
  g2 <- ggplot(pp %>% mutate(label = factor(label, levels = c("Direct (infant)",
                                                              "Total (infant)",
                                                              "Toddler",
                                                              "Adult & older adult"))), aes(coverage, est, colour=label, group=label)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    geom_line(position=position_dodge(4), linewidth=1.2) +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=5, position=position_dodge(4)) +
    geom_point(size=3, position=position_dodge(4)) +
    facet_wrap(~ panel, ncol=1, scales="free_y",
               labeller = as_labeller(c(
                 "Infants"                 = "B) Infants",
                 "Toddlers"                = "D) Toddlers",
                 "Adults & older adults"   = "F) Adults & older adults"))) +
    scale_color_nejm()+
    labs(x="Infant mAb coverage (%)", y="Protection (%) [95% CI]", colour=NULL) +
    theme_pub()
  ggsave(file.path(out_dir,"FIG2_protection.pdf"), g2, width=9, height=4.6)
  cat("[fig2] written\n")
}


icol <- c(infant="#D55E00", toddler="#009E73", adult="#0072B2",
          elderly="#CC79A7", staff="#56B4E9")   # Okabe-Ito
cov_of <- function(f) as.integer(sub(".*cov(\\d{3})\\.rds$", "\\1", basename(f)))
setting_of <- function(s) ifelse(grepl("^daycare", s), "daycare",
                                 ifelse(grepl("^household", s), "household", "community"))
have_gg <- requireNamespace("ggplot2",  quietly = TRUE); if (have_gg) library(ggplot2)
have_cl <- requireNamespace("circlize", quietly = TRUE)
files <- list.files(in_dir, "^RSV_rawsims_cov\\d{3}\\.rds$", full.names = TRUE)
if (length(files) == 0) stop(sprintf("No packs in '%s/'.", in_dir))
# ---- accumulate flows (staff excluded from household) ----------------------
flow_raw <- list(); npacks_by_cov <- integer(0)
for (cv in CONTRAST_COVS) {
  f <- files[cov_of(files) == cv]; if (length(f) == 0) next
  packs <- readRDS(f[1]); np <- length(packs); npacks_by_cov[as.character(cv)] <- np
  acc <- list()
  for (r in seq_len(np)) {
    af <- packs[[r]]$agents_final; lg <- packs[[r]]$infection_log
    if (is.null(lg) || nrow(lg) == 0) next
    rmap <- af$role; names(rmap) <- as.character(af$global_id)
    ir <- ifelse(is.na(lg$infector), "community", rmap[as.character(lg$infector)])
    ee <- lg$role_at_infection; st <- setting_of(lg$source)
    if (MERGE_STAFF_INTO_ADULT)   { ir[ir=="staff"]   <- "adult"; ee[ee=="staff"]   <- "adult" }
    if (MERGE_ELDERLY_INTO_ADULT) { ir[ir=="elderly"] <- "adult"; ee[ee=="elderly"] <- "adult" }  # <-- NEW
    ok <- ee %in% ROLES & ir %in% ROLES & st %in% SETTINGS_KEEP
    if (EXCLUDE_STAFF_IN_HH) ok <- ok & !(st=="household" & (ir=="staff" | ee=="staff"))
    if (!any(ok)) next
    key <- paste(ir[ok], ee[ok], st[ok], sep="|"); tb <- table(key)
    for (nm in names(tb)) acc[[nm]] <- (if (is.null(acc[[nm]])) 0L else acc[[nm]]) + tb[[nm]]
  }
  if (length(acc)) {
    p <- do.call(rbind, strsplit(names(acc), "|", fixed=TRUE))
    flow_raw[[as.character(cv)]] <- data.frame(coverage=cv, infector=p[,1], infectee=p[,2],
                                               setting=p[,3], count=as.integer(unlist(acc)),
                                               stringsAsFactors=FALSE)
  }
  rm(packs); gc(verbose=FALSE); cat(sprintf("cov %d: flows accumulated (%d packs)\n", cv, np))
}
flow <- do.call(rbind, flow_raw)
if (is.null(flow)) stop("No flows found.")
flow$per_centre <- flow$count / npacks_by_cov[as.character(flow$coverage)]
flow$grp_ws <- paste(flow$coverage, flow$setting, sep="|")
tot_ws <- tapply(flow$count, flow$grp_ws, sum); flow$pct_ws <- 100*flow$count/tot_ws[flow$grp_ws]
write.csv(flow[,c("coverage","setting","infector","infectee","count","per_centre","pct_ws")],
          file.path(out_dir,"FIG3_pairs.csv"), row.names=FALSE)

## ============================ FIG 3b: CHORD ===============================
if (have_cl) {
  library(circlize)
  panels <- expand.grid(setting=SETTINGS_KEEP, coverage=CONTRAST_COVS, stringsAsFactors=FALSE)
  panels <- panels[order(match(panels$setting,SETTINGS_KEEP), panels$coverage), ]
  pdf(file.path(out_dir,"FIG3b_chord.pdf"), width=10, height=10)
  par(mfrow=c(2,2), mar=c(1,1,2.2,1))
  for (i in seq_len(nrow(panels))) {
    st <- panels$setting[i]; cv <- panels$coverage[i]
    sub <- flow[flow$setting==st & flow$coverage==cv, ]
    rs <- intersect(ROLES, unique(c(sub$infector, sub$infectee)))
    mat <- matrix(0, length(rs), length(rs), dimnames=list(rs,rs))
    for (k in seq_len(nrow(sub))) mat[sub$infector[k], sub$infectee[k]] <-
      mat[sub$infector[k], sub$infectee[k]] + sub$count[k]
    tot <- sum(mat); if (tot == 0) { plot.new(); title(sprintf("%s - %d%% (no data)", st, cv)); next }
    np_cv  <- npacks_by_cov[as.character(cv)]
    out_pc <- rowSums(mat) / np_cv     # events CAUSED (as infector) per centre-season
    in_pc  <- colSums(mat) / np_cv     # events RECEIVED (as infectee) per centre-season
    circos.clear()
    circos.par(gap.degree=12, start.degree=90, points.overflow.warning=FALSE,
               canvas.xlim=c(-1.35,1.35), canvas.ylim=c(-1.35,1.35))
    chordDiagram(
      mat, grid.col=icol[rs], transparency=0.25,
      directional=1, direction.type=c("diffHeight","arrows"),
      link.arr.type="big.arrow", diffHeight=0.04,
      annotationTrack="grid", preAllocateTracks=list(track.height=0.18),
      link.sort=TRUE, link.largest.ontop=TRUE)
    # labels: on the preallocated outer track, radial, clear of the ring & arcs
    circos.trackPlotRegion(track.index=1, bg.border=NA, panel.fun=function(x,y){
      s  <- get.cell.meta.data("sector.index"); xc <- get.cell.meta.data("xcenter")
      fmt <- function(v) if (v >= 10) sprintf("%.0f", v) else sprintf("%.1f", v)
      circos.text(xc, 0.55,                       # further out on the label track
                  sprintf("%s\n%s/%s", s, fmt(out_pc[s]), fmt(in_pc[s])),
                  facing="clockwise", niceFacing=TRUE, adj=c(0.5, 0.5),
                  cex=0.9, font=2)
    })
    title(sprintf("%s  -  %d%% coverage", st, cv), cex.main=1.0, font.main=2)
  }
  circos.clear()
  par(fig=c(0,1,0,1), oma=c(1.5,0,2,0), mar=c(0,0,0,0), new=TRUE)
  plot(0,0,type="n", axes=FALSE, ann=FALSE)
  legend("top", legend=names(icol), pt.bg=icol, pch=22, pt.cex=1.9, col="grey40",
         horiz=TRUE, bty="n", title="Infector", cex=0.95)
  mtext("Who infects whom (chord, coloured by infector). Sector label = role + caused/received transmission events per centre-season. Household excludes staff.",
        side=1, outer=TRUE, cex=0.72, col="grey35")
  dev.off()
  cat("[fig3b] chord written\n")
} else cat("[fig3b] needs circlize; install.packages('circlize').\n")


if (have_gg && all(CONTRAST_COVS %in% flow$coverage)) {
  lo <- min(CONTRAST_COVS); hi <- max(CONTRAST_COVS)
  key <- function(d) paste(d$setting,d$infector,d$infectee,sep="|")
  d0 <- flow[flow$coverage==lo,]; d1 <- flow[flow$coverage==hi,]
  ks <- union(key(d0),key(d1))
  gp <- function(d,k){ v<-d$per_centre[match(k,key(d))]; v[is.na(v)]<-0; v }
  con <- data.frame(k=ks, e0=gp(d0,ks), e1=gp(d1,ks), stringsAsFactors=FALSE)
  pr <- do.call(rbind, strsplit(con$k,"|",fixed=TRUE))
  con$setting<-pr[,1]; con$infector<-pr[,2]; con$infectee<-pr[,3]
  con <- con[con$e0>=CONTRAST_MIN_EVENTS | con$e1>=CONTRAST_MIN_EVENTS, ]
  con$pct_reduction <- ifelse(con$e0>0, 100*(con$e0-con$e1)/con$e0, NA_real_)
  con$pathway <- sprintf("%s -> %s", con$infector, con$infectee)   # ASCII
  con$setting <- factor(con$setting, levels=SETTINGS_KEEP)
  con$pathway_key <- paste(con$setting, con$pathway, sep="___")
  ord_df <- unique(con[, c("pathway_key","setting","e0")])
  ord_df <- ord_df[order(ord_df$setting, -ord_df$e0), ]
  lv <- ord_df$pathway_key
  con$pathway <- factor(con$pathway_key, levels = lv, labels = sub("^.*___", "", lv))
  write.csv(con[,c("setting","infector","infectee","e0","e1","pct_reduction")],
            file.path(out_dir,"FIG3c_coverage_contrast.csv"), row.names=FALSE)

  fmt_num <- function(v) ifelse(abs(v) >= 10, sprintf("%.1f", v), sprintf("%.2f", v))
  con$disp0 <- fmt_num(con$e0)
  con$disp1 <- fmt_num(con$e1)
  con$show_change <- con$disp0 != con$disp1       # only "real" (visible) changes get an arrow
  dir_word <- ifelse(con$e1 < con$e0, "lower", "higher")
  con$lbl <- ifelse(con$show_change,
                    sprintf("%s -> %s\n(%.0f%% %s)", con$disp0, con$disp1,
                            abs(con$pct_reduction), dir_word),
                    sprintf("%s / %s", con$disp0, con$disp1))
   MIN_GAP_FRAC <- 0.035
  con$panel_max <- ave(pmax(con$e0, con$e1), con$setting, FUN = max)
  con$plot_e0 <- con$e0; con$plot_e1 <- con$e1
  ch <- con$show_change
  lo_val <- pmin(con$e0[ch], con$e1[ch]); hi_val <- pmax(con$e0[ch], con$e1[ch])
  mid     <- (lo_val + hi_val) / 2
  min_gap <- con$panel_max[ch] * MIN_GAP_FRAC
  half    <- pmax(hi_val - lo_val, min_gap) / 2
  plot_lo <- mid - half; plot_hi <- mid + half
  con$plot_e0[ch] <- ifelse(con$e0[ch] <= con$e1[ch], plot_lo, plot_hi)
  con$plot_e1[ch] <- ifelse(con$e1[ch] <  con$e0[ch], plot_lo, plot_hi)
  con_change   <- con[con$show_change, ]
  con_nochange <- con[!con$show_change, ]

  change_layer <- if (nrow(con_change) > 0)
    geom_segment(data=con_change, aes(x=pathway, xend=pathway, y=plot_e0, yend=plot_e1 * 1.005),
                 arrow=arrow(length=unit(2.2,"mm"), type="closed"),
                 colour="firebrick", linewidth=1) else NULL
  nochange_layer <- if (nrow(con_nochange) > 0)
    geom_segment(data=con_nochange, aes(x=pathway, xend=pathway, y=plot_e0, yend=plot_e1 * 1.01),
                 colour="grey75", linewidth=1) else NULL
  g <- ggplot(con) +
    change_layer +
    geom_point(aes(x=pathway, y=plot_e0, fill="0%"),   shape=21, size=2.5, colour="grey30") +
    geom_point(aes(x=pathway, y=plot_e1, fill="100%"), shape=21, size=2.5, colour="grey30") +
    geom_text(aes(x=pathway, y=pmax(plot_e0,plot_e1), label=lbl),
              vjust=-0.3, size=3.5, colour="grey15", lineheight=0.9) +
    facet_wrap(~ setting, scales="free", nrow=1,
               labeller = as_labeller(c(
                 daycare   = "A) Daycare",
                 household = "B) Household"))) +
    scale_fill_manual(values=c("0%"="white","100%"="#222222"), name="Coverage") +
    scale_y_continuous(expand=expansion(mult=c(0.04,0.22))) +
    labs(x=NULL, y="Transmission events per cohort-season") +
    theme_bw(base_size=12) +
    theme(legend.position="top",
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x=element_text(angle=40, hjust=1, size=13),
          axis.text.y=element_text(size=15),
          axis.title = element_text(size = 20),
          plot.title=element_text(face="bold"),
          strip.text = element_text(face="bold", size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          plot.caption=element_text(size=7, colour="grey45", hjust=0))
  ggsave(file.path(out_dir,"FIG3c_contrast.pdf"), g, width=18, height=8.5)
  cat("[fig3c] vertical contrast written\n")
}
cat("\nDone. FIG3b/3c + CSVs in '", out_dir, "/'.\n", sep="")


row_lv <- c("Infants","Toddlers","Adults & older adults")
cov_cols <- c(
  "0%"   = "#000000", "20%"  = "#0072B2", "40%"  = "#56B4E9",
  "60%"  = "#009E73", "80%"  = "#E69F00", "100%" = "#D55E00")
if (have_gg) {
  # ---------------------------- LEFT: prevalence -----------------------------
  panel_titles <- c(inf="Infants", tod="Toddlers", adu="Adults & older adults")
  nd <- ncol(curve_sys[[1]]$inf)
  f1 <- do.call(rbind, lapply(names(curve_sys), function(k) {
    cs <- curve_sys[[k]]
    do.call(rbind, lapply(c("inf","tod","adu"), function(s) {
      m <- cs[[s]]
      data.frame(day=seq_len(nd), coverage=as.integer(k), panel=s,
                 mean=colMeans(m), lo=apply(m,2,quantile,0.025), hi=apply(m,2,quantile,0.975),
                 row.names=NULL) })) }))
  f1$panel_lab <- factor(panel_titles[f1$panel], levels = row_lv)   # <- matches row_lv order
  cvs <- sort(unique(f1$coverage)); f1$cov_lab <- factor(paste0(f1$coverage,"%"), levels=paste0(cvs,"%"))


  y_ann <- max(f1$hi[f1$panel_lab == "Infants"], na.rm = TRUE) * 0.9
  ann_df <- data.frame(
    panel_lab = factor("Infants", levels = levels(f1$panel_lab)),
    x_text   = 93 + 45, y_text   = y_ann,
    x_arrow0 = 93 + 42, y_arrow0 = y_ann,
    x_arrow1 = 95,      y_arrow1 = y_ann
  )
  vline_df <- data.frame(panel_lab = factor("Infants", levels = levels(f1$panel_lab)))

  g1 <- ggplot(f1, aes(day, mean, colour=cov_lab, fill=cov_lab, group=cov_lab)) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.12, colour=NA) +
    geom_line(linewidth=1) +
    geom_vline(data = vline_df, aes(xintercept = 93),
               inherit.aes = FALSE, linetype = "dashed", colour = "grey30") +
    geom_segment(data = ann_df, aes(x=x_arrow0, y=y_arrow0, xend=x_arrow1, yend=y_arrow1),
                 inherit.aes = FALSE, arrow = arrow(length = unit(2,"mm"), type="closed"),
                 colour = "black", linewidth = 0.6) +
    geom_text(data = ann_df, aes(x=x_text, y=y_text, label="Catch-up dose"),
              inherit.aes = FALSE, hjust = 0, size = 3.2, colour = "black") +
    facet_wrap(~ panel_lab, ncol=1, scales="free_y",
               labeller = as_labeller(c(
                 "Infants"                 = "A) Infants",
                 "Toddlers"                = "B) Toddlers",
                 "Adults & older adults"   = "C) Adults & older adults"))) +
    scale_x_continuous(breaks=month_starts, labels=month_labels, expand=c(0.01,0)) +
    scale_colour_manual(values = cov_cols, name = "mAb coverage") +
    scale_fill_manual(values = cov_cols, name = "mAb coverage") +
    labs(x="Month", y=y_lab) +
    theme_pub() + theme(strip.text=element_text(hjust=0), legend.position = "right")

  write.csv(f1[,c("coverage","panel","day","mean","lo","hi")],
            file.path(out_dir,"FIG1_prevalence_system.csv"), row.names=FALSE)

   ARd <- function(cn,g){ t<-sum(cn[[paste0("n_",g)]]); if(t>0) sum(cn[[paste0("cd_",g)]])/t else NA_real_ }
  base_cn <- cnt_list[[as.character(min(covs))]]
  prot_point <- function(cn, base) c(
    direct=1-ARd(cn,"vaxinf")/ARd(cn,"unvaxinf"), total=1-ARd(cn,"vaxinf")/ARd(base,"unvaxinf"),
    tod=1-ARd(cn,"tod")/ARd(base,"tod"),
    adulteld=1-ARd(cn,"adulteld")/ARd(base,"adulteld"))                # <-- combined, was adult= + eld=
  boot_prot <- function(cn, base, B=B_BOOT) {
    est <- prot_point(cn, base); M <- matrix(NA_real_, B, length(est))
    for (b in seq_len(B)) {
      cb <- cn[sample.int(nrow(cn), replace=TRUE),,drop=FALSE]
      bb <- base[sample.int(nrow(base), replace=TRUE),,drop=FALSE]
      M[b,] <- prot_point(cb, bb) }
    data.frame(metric=names(est), est=est, lo=apply(M,2,quantile,0.025,na.rm=TRUE),
               hi=apply(M,2,quantile,0.975,na.rm=TRUE), row.names=NULL)
  }
  prot_df <- do.call(rbind, lapply(covs, function(cv){
    d <- boot_prot(cnt_list[[as.character(cv)]], base_cn); d$coverage <- cv; d }))
  prot_df$lo[prot_df$coverage == min(covs)] <- 0
  prot_df$hi[prot_df$coverage == min(covs)] <- 0

   lab <- c(direct="Direct (infant)", total="Total (infant)", tod="Toddler",
           adulteld="Adult & older adult")
  pnl <- c(direct="Infants", total="Infants",
           tod="Toddlers",
           adulteld="Adults & older adults")                          # <-- ONE entry now, not two
  prot_df$label <- lab[prot_df$metric]; prot_df$panel <- pnl[prot_df$metric]
  write.csv(prot_df, file.path(out_dir,"FIG2_protection.csv"), row.names=FALSE)

  pp <- prot_df[!is.na(prot_df$est), ]; pp[,c("est","lo","hi")] <- pp[,c("est","lo","hi")]*100
  pp$panel <- factor(pp$panel, levels = row_lv)   # <- SAME row order as g1
  pp$label <- factor(pp$label, levels=c("Direct (infant)","Total (infant)","Toddler","Adult & older adult"))

  g2 <- ggplot(pp, aes(coverage, est, colour=label, group=label)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey60") +
    geom_line(position=position_dodge(4), linewidth=1.2) +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=5, position=position_dodge(4)) +
    geom_point(size=3, position=position_dodge(4)) +
    facet_wrap(~ panel, ncol=1, scales="free_y",
               labeller = as_labeller(c(
                 "Infants"                 = "D) Infants",
                 "Toddlers"                = "E) Toddlers",
                 "Adults & older adults"   = "F) Adults & older adults"))) +
    scale_x_continuous(
      breaks = c(0, 20, 40, 60, 80, 100),
      labels = c(0, 20, 40, 60, 80, 100)
    ) +
    scale_color_nejm() +
    labs(x="Infant mAb coverage (%)", y="Protection (%) [95% CI]", colour=NULL) +
    theme_pub() + theme(legend.position="right")

   if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    g_combined <- g1 | g2

    ggsave(file.path(out_dir,"FIG1_2_combined.pdf"), g_combined, width=15, height=8)
    cat("[fig1+2] combined figure written\n")
  } else {
    ggsave(file.path(out_dir,"FIG1_prevalence_system.pdf"), g1, width=8.5, height=7)
    ggsave(file.path(out_dir,"FIG2_protection.pdf"), g2, width=8.5, height=7)
    cat("[fig1] and [fig2] written separately -- install.packages('patchwork') to combine them.\n")
  }
}


#source("LRTI_Hospitalization.R")
