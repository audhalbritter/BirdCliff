# traits_raw %>% distinct(Taxon)
#
# # choosing species and traits
# traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count() %>%
#   filter(n > 3)
#
# bryo_traits_raw %>%
#   distinct(Gradient, Site, Taxon) %>%
#   group_by(Gradient, Taxon) %>%
#   count()


combine_traits <- function(traits_raw, bryo_traits_raw){

  # cobmine species vascular and bryophytes
  all_traits <- bind_rows(traits_raw %>%
                            filter(Taxon %in% c("luzula confusa", "salix polaris", "cerastium arcticum"),
                                   trait_trans %in% c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil")),
                          bryo_traits_raw %>%
                            filter(Taxon %in% c("aulacomnium turgidum", "hylocomium splendens", "sanionia uncinata"),
                                   trait_trans %in% c("Shoot_Length_cm_log", "Shoot_ratio", "SSL_cm_g", "WHC_g_g", "P_percent"))) %>%
    mutate(trait_trans = factor(trait_trans, levels = c("Plant_Height_cm_log", "Leaf_Area_cm2_log", "LDMC", "N_percent", "dN15_permil", "Shoot_Length_cm_log", "Shoot_ratio", "SSL_cm_g", "WHC_g_g", "P_percent")),
           Taxon = factor(Taxon, levels = c("cerastium arcticum", "luzula confusa", "salix polaris", "aulacomnium turgidum", "hylocomium splendens", "sanionia uncinata"))) %>%
    mutate(GS = paste0(Gradient, Site))

  return(all_traits)
}

make_ind_sp_plot <- function(all_traits){

  v <- fancy_trait_name_dictionary(all_traits) %>%
    filter(Functional_group == "vascular") %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))


  b <- fancy_trait_name_dictionary(all_traits) %>%
    filter(Functional_group == "bryophyte") %>%
    ggplot(aes(x = Elevation_m, y = value_trans, colour = Gradient)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_colour_manual(values = c("green4", "grey"), labels = c("Bird cliff", "Reference")) +
    facet_grid(trait_fancy ~ Taxon, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(aspect.ratio = 0.8,
          plot.margin = margin(0.5, 0.5, 0.5, 0.5))

  ind_sp_traits <- v + b + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  return(ind_sp_traits)

}



#cerastium arcticum 4 Height, Area, N, N15
salix polaris 1 Height
luzula confusa

#sanionia uncinata   4 P shoot length shoot ratio ssl
#aulacomnium turgidum 4 shoot length shoot ratio ssl shc
hylocomium splendens 3 P shoot ratio WHC

1 Plant_Height_cm_log
2 Leaf_Area_cm2_log
3 LDMC
4 N_percent
5 dN15_permil
6 P_percent
7 Shoot_Length_cm_log
8 Shoot_ratio
9 SSL_cm_g
10 WHC_g_g

dd <- all_traits %>%
  filter(trait_trans == "WHC_g_g",
         Taxon == "hylocomium splendens")
fit <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), na.action = "na.fail", REML = FALSE, data = dd)
model.set = dredge(fit, rank = "AICc", extra = "R^2")
model.set

sal.pol A: G+E, LDMC: G+E, N: G, N15: G
#cer.arc LDMC: Null
luz.conf H: G+E, A: G+E, LDMC: GxE, N: G, N15: G

san.unc WHC: G
aul.tur P: Null
hyl spl SL: Null, SSL: Null


all_traits %>% distinct(Taxon)

# run models for single species
all_traits %>%
  filter(Taxon %in% c("salix polaris", "luzula confusa")) %>%
  group_by(Taxon, trait_trans) %>%
  nest(data = -c(Taxon, trait_trans)) %>%
  mutate(estimate = map(data, ~{
    mod <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), data =  .x)
    estimates = broom.mixed::tidy(mod)
  })) %>%
  unnest(estimate) %>%
  mutate(r = map(data, ~{
    mod <- lmer(value_trans ~ Gradient * Elevation_m + (1|GS), data =  .x)
    r = as.numeric(r.squaredGLMM(mod))
  })) %>%
  unnest_wider(col = r) %>%
  select(trait_trans, "Rm" = "...1", "Rc" = "...2")





all_traits %>%
  filter(Taxon == "salix polaris") %>%
  group_by(trait_trans) %>%
  nest(data = -c(trait_trans)) %>%
  mutate(model.set = map(data, ~{
    mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
    model.set = dredge(mod, rank = "AICc", extra = "R^2")
    model.set
  })) %>%
  unnest(model.set)





all_traits %>%
  filter(Taxon == "luzula confusa") %>%
  group_by(trait_trans) %>%
  nest(data = -c(trait_trans)) %>%
  mutate(mod = map(data, ~{
    mod <- lmer(mean ~  Gradient * Elevation_m + (1|GS), data = .x)
  }))


all_traits %>%
  filter(Taxon == "luzula confusa") %>%
  group_by(trait_trans) %>%
  nest(data = -c(trait_trans)) %>%
  mutate(estimate = map(data, ~{
    mod <- lmer(mean ~ Gradient + (1|Site), data = .x)
    estimates = broom.mixed::tidy(mod)
  })) %>%
  unnest(estimate)

