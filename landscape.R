# ..
require(ggplot2)
require(scales)
require(data.table)

d <- fread('data/Pairwise_cultivations_100x_270319.csv')
d <- d[, .(Co_fermentations, Yield_Average, Yield_SD)]
names(d) <- c('comp', 'yld', 'yld.se')


# compute the function as log[Yield_(yeast + community) / Yield_(yeast)]
d[, F := log(yld / d[comp =='yeast']$yld)]
d[, F.sd := sqrt((yld.se/yld)^2 + (d[comp =='yeast']$yld.se/d[comp =='yeast']$yld)^2)]
d.yeast = d[comp == 'yeast']
d = d[comp != 'yeast']
order.strain.names <- function(x) return(paste(sort(unlist(strsplit(x, '-'))), collapse='-'))
d[, comp := order.strain.names(comp), by = seq_len(nrow(d))]
count.strains <- function(x) return(length(unlist(strsplit(x, '-'))))
d[, c.len := count.strains(comp), by = seq_len(nrow(d))]

# discard data for which error is too high
d = d[comp != 'sp2st2-sp3-sp4']

# discard species that have been only cocultured with yeast
d = d[comp %nin% c('sp7', 'sp2st3', 'sp8')]

strains <- unique(unlist(strsplit(d$comp, '-')))


# =============================================================================
# .. the "epistasis" map
# =============================================================================

# .. first get all backgrounds (i.e. combinations of <=5 strains, including 0)
bg <- unique(d$comp)
bg <- bg[sapply(strsplit(bg, '-'), length)<=8]
bg <- c(bg, NA)

# .. and the pairs
pair <- unique(d$comp)
pair <- pair[sapply(strsplit(pair, '-'), length)==2]

# .. get map between pairs and all combinations of the rest of strains
map <- data.table(expand.grid(bg, pair, stringsAsFactors = FALSE))
names(map) <- c('bg', 'pair')
auxfun <- function(x, y) any(unlist(strsplit(x,'-')) %in% unlist(strsplit(y,'-')))
map[, all.in := auxfun(pair, bg), by = 1:nrow(map)]
map <- map[all.in==FALSE]
map[, all.in :=NULL]

# .. background function
map[, bg.obs := d[match(bg, comp)]$F]
map[, bg.obs.se := d[match(bg, comp)]$F.sd]

map[, bg.obs := ifelse(is.na(bg), 0, bg.obs)]
map[, bg.obs.se := ifelse(is.na(bg), 0, bg.obs.se)]

# .. compute expected: observed F for BG + expected from singles
r <- map[,.(paste(sort(c(unlist(strsplit(bg, '-')),
                         unlist(strsplit(pair, '-'))[1])), collapse = '-')),
         by = 1:nrow(map)]$V1

map[, single.1 := d[match(r, comp)]$F - bg.obs]
map[, s1.e := sqrt(d[match(r, comp)]$F.sd^2 + bg.obs.se^2)]

r <- map[,.(paste(sort(c(unlist(strsplit(bg, '-')),
                         unlist(strsplit(pair, '-'))[2])), collapse = '-')),
         by = 1:nrow(map)]$V1
map[, single.2 := d[match(r, comp)]$F - bg.obs]
map[, s2.e := sqrt(d[match(r, comp)]$F.sd^2 + bg.obs.se^2)]

map[, expected.F := single.1 + single.2]
map[, expected.se := sqrt(s1.e^2 + s2.e^2)]

# .. add the observed one to the table
r <- map[,.(paste(sort(c(unlist(strsplit(bg, '-')),
                         unlist(strsplit(pair, '-')))), collapse = '-')),
         by = 1:nrow(map)]$V1

map[, observed.F := d[match(r, comp)]$F - bg.obs]
map[, observed.se := sqrt(d[match(r, comp)]$F.sd^2 + bg.obs.se^2)]

# .. compute epsilon
map[, epsilon := observed.F - expected.F]
map[, epsilon.se := sqrt(observed.se^2 + expected.se^2)]

# .. add column for epsilon^0
aux <- map[is.na(bg)]
aux <- aux[, c('pair', 'epsilon', 'epsilon.se')]
names(aux)[2:3] <- c('eps.zero', 'eps.zero.se')
map <- merge(map, aux, by = 'pair', all=TRUE)

# .. add background community size
map[, bg.size := sapply(strsplit(bg, '-'), length)]
map[, bg.size := ifelse(is.na(bg), 0, bg.size)]

# .. delta_epsilon_ABC
map[, d.epsilon.ABC := epsilon-eps.zero]
map[, d.epsilon.ABC.se := sqrt(epsilon.se^2 + eps.zero.se^2)]



# .. pairwise epistasis heatmap
map.1 <- copy(map)
map.1[, p.1 := sapply(strsplit(pair, '-'), '[',1)]
map.1[, p.2 := sapply(strsplit(pair, '-'), '[',2)]

map.2 <- copy(map)
map.2[, p.1 := sapply(strsplit(pair, '-'), '[',2)]
map.2[, p.2 := sapply(strsplit(pair, '-'), '[',1)]

aux <- rbind(map.1, map.2)

k <- subset(aux, is.na(bg))
k <- k[,.(sum(abs(epsilon))), by = p.1]
k <- k[order(V1, decreasing=TRUE)]
k <- k$p.1

aux[, p.1 := ordered(p.1, rev(k))]
aux[, p.2 := ordered(p.2, rev(k))]

clr <- brewer_pal('div', 5)(9)

k <- subset(aux, is.na(bg) & p.1<p.2)
k[, eps.signif := ifelse(epsilon-(2*epsilon.se)<0 & epsilon+(2*epsilon.se)>0,
                         '', '*')]

p <- 'results/pairwise_epistasis_checkerbox.pdf'
pdf(p, 4.5, 4)
ggplot(k, aes(x=p.1, y=p.2, fill = epsilon)) +
  geom_tile(colour = NA) +
  geom_text(aes(label = eps.signif), size=10) +
  scale_fill_gradient2(low = clr[1], mid = clr[5], high=clr[9], midpoint=0) +
  theme_bw() + labs(title = 'Pairwise epistasis', x = 'Species 1', y = 'Species 2') +
  theme(panel.grid = element_blank())
dev.off()



# =============================================================================
# .. the "structure-function landscape"
#  =============================================================================

# .. first build all possible trajectories (structure-function landscape)
strfun <- data.table(t(combn(d$comp, 2)))
strfun[, len.v1 := d[match(strfun$V1, d$comp), c.len]]
strfun[, len.v2 := d[match(strfun$V2, d$comp), c.len]]
strfun.V1 <- with(strfun, ifelse(len.v1<len.v2, V1, V2))
strfun.V2 <- with(strfun, ifelse(len.v1<len.v2, V2, V1))
strfun$V1 <- strfun.V1
strfun$V2 <- strfun.V2
strfun <- strfun[,1:2]
strfun[, len.v1 := d[match(strfun$V1, d$comp), c.len]]
strfun[, len.v2 := d[match(strfun$V2, d$comp), c.len]]
strfun <- strfun[len.v2 == (len.v1+1)]

auxfun <- function(x, y) all(unlist(strsplit(x,'-')) %in% unlist(strsplit(y,'-')))
strfun[, all.in := auxfun(V1, V2), by = 1:nrow(strfun)]
strfun <- strfun[all.in==TRUE]

# .. get real Vj (for actual landscape) and standard errors
strfun$F.real.V1 <- d$F[match(strfun$V1, d$comp)]
strfun$F.real.V2 <- d$F[match(strfun$V2, d$comp)]

strfun$se.real.V1 <- d$F.sd[match(strfun$V1, d$comp)]
strfun$se.real.V2 <- d$F.sd[match(strfun$V2, d$comp)]

strfun <- unique(strfun)


# .. get str-fun landscape if everything were additive **
aux <- d[c.len==1]

auxfun <- function(x) sum(aux$F[match(unlist(strsplit(x, '-')), aux$comp)])
strfun[, F.1.V1 := auxfun(V1), by = 1:nrow(strfun)]
strfun[, F.1.V2 := auxfun(V2), by = 1:nrow(strfun)]

# errors **
auxfun <- function(x) sqrt(sum(aux$F.sd[match(unlist(strsplit(x, '-')), aux$comp)]^2))
strfun[, F.1.V1.se := auxfun(V1), by = 1:nrow(strfun)]
strfun[, F.1.V2.se := auxfun(V2), by = 1:nrow(strfun)]


# .. get str-fun landscape if truncating at pairwise
aux <- map[is.na(bg)]
pairwise <- function(x) {
    r <- unlist(strsplit(x, '-'))
    if(length(r)>1) {
        r <- apply(t(combn(r, 2)), 1, paste, collapse = '-')
        r.val <- sum(aux$epsilon[match(r, aux$pair)], na.rm=TRUE) ## here na.rm is wrong actually, so it will give us only an estimate based on known epsilon values
        r.err  <- sqrt(sum(aux$epsilon.se[match(r, aux$pair)]^2, na.rm=TRUE))
    } else {
        r.val <- 0
        r.err  <-  0
    }
    return(list(r.val, r.err))
}

strfun[, c('pair.int', 'pair.se'):= pairwise(V1), by = 1:nrow(strfun)]
strfun[, F.2.V1 := F.1.V1 + pair.int]
strfun[, F.2.V1.se := sqrt(F.1.V1.se^2 + pair.se^2)]
strfun[, c('pair.int', 'pair.se'):= NULL]

strfun[, c('pair.int', 'pair.se'):= pairwise(V2), by = 1:nrow(strfun)]
strfun[, F.2.V2 := F.1.V2 + pair.int]
strfun[, F.2.V2.se := sqrt(F.1.V2.se^2 + pair.se^2)]
strfun[, c('pair.int', 'pair.se'):= NULL]

# strfun[, F.2.V1 := F.1.V1 + pairwise.int(V1), by = 1:nrow(strfun)]
# strfun[, F.2.V2 := F.1.V2 + pairwise.int(V2), by = 1:nrow(strfun)]


# .. get Vj taking three-wise interactions into account
# .. first compute third-order interactions
aux <- strfun[len.v2 == 3][,c('V2', 'F.2.V2', 'F.2.V2.se', 'F.real.V2', 'se.real.V2')]
aux <- unique(aux)
aux[, eps.3 := F.real.V2 - F.2.V2]
aux[, eps.3.se := sqrt(se.real.V2^2 + F.2.V2.se^2)]

auxfun <- function(x) {
    r <- unlist(strsplit(x, '-'))
    if(length(r)>2) {
        r <- apply(t(combn(r, 3)), 1, paste, collapse = '-')
        r.val <- sum(aux$eps.3[match(r, aux$V2)], na.rm=TRUE)
        r.err  <- sqrt(sum(aux$eps.3.se[match(r, aux$V2)]^2, na.rm=TRUE))
    } else {
        r.val <- 0
        r.err  <-  0
    }
    return(list(r.val, r.err))
}


strfun[, c('pair.int', 'pair.se'):= auxfun(V1), by = 1:nrow(strfun)]
strfun[, F.3.V1 := F.2.V1 + pair.int]
strfun[, F.3.V1.se := sqrt(F.2.V1.se^2 + pair.se^2)]
strfun[, c('pair.int', 'pair.se'):= NULL]

strfun[, c('pair.int', 'pair.se'):= auxfun(V2), by = 1:nrow(strfun)]
strfun[, F.3.V2 := F.2.V2 + pair.int]
strfun[, F.3.V2.se := sqrt(F.2.V2.se^2 + pair.se^2)]
strfun[, c('pair.int', 'pair.se'):= NULL]


# .. add singles to the str-fun landscape
aux <- strfun[len.v1==1]
aux <- aux[,c(2,1, 4,3,5, 7,6, 9,8, 11,10, 13,12, 16,17, 14,15, 20, 21, 18, 19)]
names(aux) <- names(strfun)
aux$V1 <- '-'
aux$len.v1 <- 0
aux$F.real.V1 <- 0
aux$se.real.V1 <- 0
aux$F.1.V1 <- 0
aux$F.1.V1.se <- 0
aux$F.2.V1 <- 0
aux$F.3.V1 <- 0

aux <- unique(aux)

strfun <- rbind(strfun, aux)


# .. variation upon species removal
strfun[, vari := F.real.V2 - F.real.V1]


# =============================================================================
# .. plots
# =============================================================================

p <- 'results/fig_2A.pdf'
pdf(p, 4, 4, family = 'CM Sans')
ggplot(d, aes(x = c.len, y = F)) +
geom_point(data = strfun, colour = 'blue', alpha=.4, size = 1,
              aes(x = len.v2, y = F.1.V2)) +
geom_segment(data = strfun, colour = 'gray30', alpha=.4,
                mapping = aes(x = len.v1, xend = len.v2,
                              y = F.1.V1, yend = F.1.V2)) +
labs(x = 'Number of species', y='F') + ylim(-1.1, .1) +
theme_bw() + labs(title = 'Order 1 (only single effects)') +
theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)


p <- 'results/fig_expected_pairwise.pdf'
pdf(p, 4, 4, family = 'CM Sans')
ggplot(d, aes(x = c.len, y = F)) +
geom_point(data = strfun, colour = 'gray30', alpha=.4, size = 1,
              aes(x = len.v2, y = F.2.V2)) +
geom_segment(data = strfun, colour = 'red', alpha=.4,
                mapping = aes(x = len.v1, xend = len.v2,
                              y = F.2.V1, yend = F.2.V2)) +
labs(x = 'Number of species', y='Vj') + ylim(-1.1, .1) +
theme_bw() + labs(title = 'Order 2 (including pairwise effects)') +
theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)


k <- strfun[, .(mean.bylen = mean(F.real.V2),
                se.bylen = sd(F.real.V2)),
            by=len.v2]

p <- 'results/fig_2B.pdf'
pdf(p, 4, 4)
ggplot(d, aes(x = c.len, y = F)) +
  geom_point(data = strfun, colour = 'gray30', alpha=.4, size = 1,
             aes(x = len.v2, y = F.real.V2)) +
  geom_segment(data = strfun, colour = 'gray30', alpha=.5,
               mapping = aes(x = len.v1, xend = len.v2,
                             y = F.real.V1, yend = F.real.V2)) +
  geom_point(data=k, aes(len.v2, mean.bylen), color='red2', size=2) +
  geom_line(data=k, aes(len.v2, mean.bylen), color='red2', size=1) +
  labs(x = 'Number of species', y='Vj') + ylim(-1.1, .1) +
  theme_bw() + labs(title = 'Observed landscape') +
  theme(panel.grid = element_blank())
dev.off()

p <- 'results/fig_3A_subset_with_sp1.pdf'
pdf(p, 4, 4)
ggplot(d, aes(x = c.len, y = F)) +
    geom_point(data = strfun[V1 %like% 'sp1' & V2 %like% 'sp1'],
               colour = 'gray30', alpha=.4, size = 1,
               aes(x = len.v2, y = F.real.V2)) +
    geom_segment(data = strfun[V1 %like% 'sp1' & V2 %like% 'sp1'],
                 colour = 'gray30', alpha=.5,
                 mapping = aes(x = len.v1, xend = len.v2,
                               y = F.real.V1, yend = F.real.V2)) +
    labs(x = 'Number of species', y='Vj') + ylim(-1.1, .1) +
    theme_bw() + labs(title = 'Observed landscape') +
    theme(panel.grid = element_blank())
dev.off()

p <- 'results/fig_3A_subset_without_sp1.pdf'
pdf(p, 4, 4)
ggplot(d, aes(x = c.len, y = F)) +
    geom_point(data = strfun[!(V1 %like% 'sp1') & !(V2 %like% 'sp1')],
               colour = 'gray30', alpha=.4, size = 1,
               aes(x = len.v2, y = F.real.V2)) +
    geom_segment(data = strfun[!(V1 %like% 'sp1') & !(V2 %like% 'sp1')],
                 colour = 'gray30', alpha=.5,
                 mapping = aes(x = len.v1, xend = len.v2,
                               y = F.real.V1, yend = F.real.V2)) +
    labs(x = 'Number of species', y='Vj') + ylim(-1.1, .1) +
    theme_bw() + labs(title = 'Observed landscape') +
    theme(panel.grid = element_blank())
dev.off()


p <- 'results/fig_2A_inset.pdf'
pdf(p, 3, 3)
ggplot(k, aes(x = len.v2, y = mean.bylen)) +
  geom_pointrange(aes(ymin=mean.bylen-se.bylen, ymax=mean.bylen+se.bylen)) +
  geom_hline(yintercept=0, col='gray', linetype=2) +
  labs(x = 'Number of species', y='F') +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()




# .. epistasis in pairs as a function of single background strains (3-body)
aux <- map[bg %in% strains | is.na(bg)]
aux[, bg := ifelse(is.na(bg), '-', bg)]
aux[, bg := ordered(bg, c('-', sort(strains)))]

p <- 'results/fig_3B.pdf'
pdf(p, 6, 3, family = 'CM Sans')
ggplot(subset(aux, bg!='-'), aes(as.factor(bg), d.epsilon.ABC, color=pair)) +
    geom_hline(yintercept = 0, col='gray') +
    geom_pointrange(aes(ymin=d.epsilon.ABC - d.epsilon.ABC.se,
                        ymax=d.epsilon.ABC + d.epsilon.ABC.se),
                    position=position_jitter(width=0.15)) +
    theme_bw() + labs(title = 'Pairs cocultured with one background strain',
                      x = 'Background strain') +
    theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)


# .. strength of interactions in 3-body
k <- cbind(do.call(rbind, strsplit(aux$pair, '-')), as.character(aux$bg))
k <- t(apply(k, 1, sort))
k <- apply(k, 1, paste, collapse = '-')
k <- strfun[match(k, V2)]

aux <- cbind(aux, k[, c('F.1.V2', 'F.1.V2.se',
                        'F.2.V2', 'F.2.V2.se',
                        'F.real.V2', 'se.real.V2')])
aux[, eps.2 := (F.2.V2 - F.1.V2)]
aux[, eps.2.se := sqrt(F.2.V2.se^2 + F.1.V2.se^2)]


p <- 'results/fig_3C.pdf'
pdf(p, 5, 4, family = 'CM Sans')
ggplot(aux, aes(eps.2, d.epsilon.ABC)) +
geom_errorbar(aes(ymin=d.epsilon.ABC - d.epsilon.ABC.se,
                  ymax=d.epsilon.ABC + d.epsilon.ABC.se),
              size=.5, width = .01, colour = 'gray') +
geom_errorbarh(aes(xmin=eps.2 - eps.2.se,
                   xmax=eps.2 + eps.2.se),
               size=.5, width = .01, colour = 'gray') +
geom_hline(yintercept = 0, col='gray', linetype = 2) +
geom_vline(xintercept = 0, col='gray', linetype = 2) +
geom_point(size = 3) +
theme_bw() + labs(title = 'Triads',
                  y = 'Higher order epistasis' ,
                  x = expression(paste(Sigma, epsilon[ij]^0))) +
theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)


p <- 'results/fig_S2.pdf'
pdf(p, 5, 4, family = 'CM Sans')
ggplot(aux, aes(F.1.V2, eps.2)) +
geom_hline(yintercept = 0, col='gray', linetype = 2) +
geom_vline(xintercept = 0, col='gray', linetype = 2) +
geom_errorbar(aes(ymin=eps.2 - eps.2.se,
                  ymax=eps.2 + eps.2.se),
              size=.5, width = 0, colour = 'gray') +
geom_errorbarh(aes(xmin=F.1.V2 - F.1.V2.se,
                   xmax=F.1.V2 + F.1.V2.se),
               size=.5, width = 0, colour = 'gray') +
geom_point(size = 3) +
scale_colour_manual(values = c('skyblue3', 'chocolate2')) +
theme_bw() + labs(title = 'Triads',
                  y = 'Pairwise interaction' ,
                  x = 'Additive expectation (1st order)') +
theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)


# .. higher order deviations as a function of community size

strfun[, eps.h := F.real.V2 - F.2.V2]
strfun[, eps.h.se := sqrt(se.real.V2^2 + F.2.V2.se^2)]
strfun[, eps.h.se := ifelse(len.v2 == 2, 0, eps.h.se)]



p <- 'results/fig_S13A.pdf'
pdf(p, 6, 4, family = 'CM Sans')
ggplot(strfun, aes(F.2.V2, F.real.V2, colour = len.v2)) +
geom_errorbarh(aes(xmin=F.2.V2 - F.2.V2.se,
                   xmax=F.2.V2 + F.2.V2.se),
               size=.5, width = 0, colour = 'gray') +
geom_errorbar(aes(ymin=F.real.V2 - se.real.V2,
                  ymax=F.real.V2 + se.real.V2),
              size=.5, width = 0, colour = 'gray') +
geom_point(size = 3) + #xlim(0, 100) + ylim(0, 50) +
theme_bw() + labs(y = 'Measured' ,
                  x = 'Predicted using pairwise model',
                  colour = 'Species in\n community') +
geom_abline(intercept = 0, slope = 1) +
theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)



p <- 'results/fig_2F.pdf'
pdf(p, 5, 4, family = 'CM Sans')
ggplot(strfun, aes(x = F.1.V2, y = F.real.V2,
                   shape = as.factor(len.v2))) +
    geom_errorbar(aes(ymin=F.real.V2 - se.real.V2,
                      ymax=F.real.V2 + se.real.V2),
                  size=.3, width = 0, colour = 'gray') +
    geom_errorbarh(aes(xmin=F.1.V2 - F.1.V2.se,
                       xmax=F.1.V2 + F.1.V2.se),
                   size=.3, height = 0, colour = 'gray') +
    geom_point(size = 3, color = 'turquoise4') +
    theme_bw() + labs(y = 'Measured' ,
                      x = 'Predicted from additive model ',
                      colour = 'Species in\n community') +
    scale_shape_manual(values = c(c(0:2, 15:18))) +
    geom_abline(intercept = 0, slope = 1) +
    theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)

aux <- strfun[, .(V2, F.real.V2, se.real.V2, F.1.V2, F.1.V2.se)]
aux <- unique(aux)

aux[, sum((F.real.V2<(F.1.V2+F.1.V2.se) & F.real.V2>(F.1.V2-F.1.V2.se)) |
          (F.1.V2<(F.real.V2+se.real.V2) & F.1.V2>(F.real.V2-se.real.V2)))]




aux <- strfun
aux[, pred.PA := F.real.V2-F.1.V2]
aux[, pred.PA.se := sqrt(se.real.V2^2 + F.1.V2.se^2)]

## aux <- aux[, c('V2', 'len.v2', 'pred.PA', 'pred.pairwise')]
## aux <- unique(aux)
## aux <- melt(aux, id.vars=c('V2', 'len.v2'), measure.vars=c('pred.PA', 'pred.pairwise'))


aux <- strfun[, c('V2', 'len.v2', 'pred.PA', 'pred.PA.se')]
aux <- unique(aux)

p <- 'results/fig_predictions_simple.pdf'
pdf(p, 5, 4)
ggplot(aux, aes(factor(len.v2), pred.PA)) +
  geom_hline(yintercept=0, color='gray', linetype=2) +
  geom_pointrange(aes(ymin=pred.PA-pred.PA.se, ymax=pred.PA+pred.PA.se), size=.2,
                  position = position_jitter(width = 0.1)) +
  ylab('|Measured - Predicted|') +
  xlab('Community size') +
  theme(panel.grid=element_blank())
dev.off()




aux <- strfun[, c('V2', 'len.v2', 'eps.h', 'eps.h.se')]
aux <- unique(aux)

p <- 'results/fig_predictions_pairwise.pdf'
pdf(p, 5, 4, family = 'CM Sans')
ggplot(subset(aux, len.v2>1), aes(factor(len.v2), eps.h)) +
  geom_hline(yintercept=0, color='gray', linetype=2) +
  geom_pointrange(aes(ymin=eps.h-eps.h.se, ymax=eps.h+eps.h.se), size=.2,
                  position = position_jitter(width = 0.1)) +
  ylab('|Measured - Predicted|') +
  xlab('Community size') +
  theme(panel.grid=element_blank())
dev.off()
embed_fonts(p, outfile = p)


# how many communities are well predicted
aux[, F.diff := F.real.V2-F.1.V2]
aux[, F.diff.se := sqrt((se.real.V2^2)+(F.1.V2.se^2))]

aux[, signif := (F.diff-F.diff.se)<0 & (F.diff+F.diff.se)>=0]

# .. interactions of higher order than 3
strfun[, eps.h3 := F.real.V2 - F.3.V2]
strfun[, eps.h3.se := sqrt(se.real.V2^2 + F.3.V2.se^2)]
strfun[, eps.h3.se := ifelse(len.v2 <= 3, 0, eps.h3.se)]

p <- 'results/fig_S13C.pdf'
pdf(p, 5, 4, family = 'CM Sans')
ggplot(strfun, aes(len.v2, abs(eps.h3), colour = V2 %like% 'P')) +
    geom_hline(yintercept = 0, col='gray', linetype = 2) +
    geom_errorbar(aes(ymin=abs(eps.h3) - eps.h3.se,
                      ymax=abs(eps.h3) + eps.h3.se),
                  size=.5, width = 0, colour = 'gray') +
    geom_point(size = 3) +
    scale_colour_manual(values = c('skyblue3', 'chocolate2')) +
    theme_bw() + labs(title = 'Communities',
                      y = 'Higher order epistasis' ,
                      x = 'Community size',
                      colour = 'P. polymyxa \nin community') +
    theme(panel.grid = element_blank())
dev.off()
embed_fonts(p, outfile = p)

