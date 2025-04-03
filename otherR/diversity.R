#test diversity along gradients
tar_target(
  name = diversity_analysis,
  command = {

    best <- diversity_grad %>%
      ungroup() %>%
      distinct(DiversityIndex) %>%
      mutate(best_model = c("G", "G", "Null", "Null"))

    # GRADIENT MODEL only richness and diverity
    nest <- diversity_grad %>%
      mutate(GS = paste0(Gradient, Site)) %>%
      filter(DiversityIndex %in% c("Richness", "Diversity")) %>%
      group_by(DiversityIndex) %>%
      nest(data = -c(DiversityIndex))

    r_square <- nest %>%
      mutate(r = map(data, ~{
        mod <- lmer(Value ~ Gradient + (1|GS), data = .x)
        r = as.numeric(r.squaredGLMM(mod))
      })) %>%
      unnest_wider(col = r) %>%
      select(DiversityIndex, "Rm" = "...1", "Rc" = "...2")

    estimate <- nest %>%
      mutate(mod = map(data, ~lmer(Value ~ Gradient + (1|GS), data = .x)),
             result = map(mod, tidy)) %>%
      unnest(result)

    # NULL MODEL evenness and sumAbundance
    nest_2 <- diversity_grad %>%
      mutate(GS = paste0(Gradient, Site)) %>%
      filter(DiversityIndex %in% c("Evenness", "sumAbundance")) %>%
      group_by(DiversityIndex) %>%
      nest(data = -c(DiversityIndex))

    r_square_2 <- nest_2 %>%
      mutate(r = map(data, ~{
        mod <- lmer(Value ~ 1 + (1|GS), data = .x)
        r = as.numeric(r.squaredGLMM(mod))
      })) %>%
      unnest_wider(col = r) %>%
      select(DiversityIndex, "Rm" = "...1", "Rc" = "...2")

    estimate_2 <- nest_2 %>%
      mutate(mod = map(data, ~lmer(Value ~ 1 + (1|GS), data = .x)),
             result = map(mod, tidy)) %>%
      unnest(result)

    diversity_output <- bind_rows(estimate, estimate_2) %>%
      select(-data, -mod) %>%
      filter(effect == "fixed") %>%
      left_join(best, by = "DiversityIndex") %>%
      left_join(bind_rows(r_square, r_square_2), by = "DiversityIndex") %>%
      select(Index = DiversityIndex, "Best model" = best_model, term, Estimate = estimate, "Standard error" = std.error, "t-value" = statistic, "Marginal R2" = Rm, "Conditional R2" = Rc)

    return(diversity_output)

  })

# test best diversity model
tar_target(
  name = div_best_model,
  command = {
    diversity_grad %>%
      mutate(GS = paste0(Gradient, Site)) %>%
      group_by(DiversityIndex) %>%
      nest(data = -c(DiversityIndex)) %>%
      mutate(model.set = map(data, ~{
        mod <- lmer(Value ~  Gradient * Elevation_m + (1|GS), REML = FALSE, na.action = "na.fail", data = .x)
        model.set = dredge(mod, rank = "AICc", extra = "R^2")
      })) %>%
      unnest(model.set)
  })
