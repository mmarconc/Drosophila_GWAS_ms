# ============================================================
# Paired Kaplan–Meier + CoxME Mixed Model (dish as random effect)
# ============================================================

library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(broom)
library(coxme)

# ------------------------------------------------------------
# 1. Load and prepare data
# ------------------------------------------------------------
df <- read_excel("Axenic.xlsx")

df <- df %>%
  mutate(
    conc = as.numeric(conc),
    dish = as.factor(dish),
    name = str_trim(as.character(name))
  ) 

df$SurvObj <- with(df, Surv(time, death))

# ------------------------------------------------------------
# 2. Summarize survival probability at time = 60
# ------------------------------------------------------------
target_time <- 60

summary_df <- df %>%
  group_by(name, conc) %>%
  group_modify(~{
    fit <- tryCatch(survfit(Surv(time, death) ~ 1, data = .x), error = function(e) NULL)
    if (is.null(fit) || length(fit$surv) == 0) {
      return(tibble(surv = 0, lower = 0, upper = 0))
    }
    s <- summary(fit, times = target_time)
    if (length(s$surv) == 0) {
      tibble(surv = 0, lower = 0, upper = 0)
    } else {
      tibble(surv = s$surv, lower = s$lower, upper = s$upper)
    }
  }) %>%
  ungroup()

# ------------------------------------------------------------
# 3. Identify base vs variant pairs (e.g., "CS" ↔ "CSA")
# ------------------------------------------------------------
summary_df <- summary_df %>%
  mutate(
    base = str_replace(name, "A$", ""),  # remove trailing "A"
    is_variant = str_detect(name, "A$")
  )

paired_bases <- summary_df %>%
  group_by(base) %>%
  filter(any(is_variant) & any(!is_variant)) %>%
  pull(base) %>%
  unique()

summary_df <- summary_df %>% filter(base %in% paired_bases)

# ------------------------------------------------------------
# 4. Mixed-effects Cox model (dish as random effect) per pair
# ------------------------------------------------------------
coxme_results <- df %>%
  mutate(
    base = str_replace(name, "A$", ""),
    is_variant = str_detect(name, "A$")
  ) %>%
  filter(base %in% paired_bases) %>%
  group_by(base) %>%
  group_modify(~{
    if (length(unique(.x$is_variant)) < 2) {
      return(tibble(p.value = NA))
    }
    
    fit <- try(coxme(Surv(time, death) ~ is_variant + conc + (1 | dish), data = .x), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      return(tibble(p.value = NA))
    }
    
    beta <- fixef(fit)["is_variantTRUE"]
    se <- sqrt(diag(vcov(fit)))["is_variantTRUE"]
    z <- beta / se
    p <- 2 * pnorm(-abs(z))
    
    tibble(p.value = p)
  }) %>%
  ungroup()

coxme_results

# ------------------------------------------------------------
# 5. Plot paired dose–response with CI and CoxME results
# ------------------------------------------------------------

summary_df$conc_f <- factor(summary_df$conc, levels = c(2, 3, 4, 5, 7, 10, 15))


p1 <- ggplot(summary_df, aes(x = conc_f, y = surv, color = name, group = name)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = name),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  # facets now follow the factor order set above
  facet_wrap(~ base, scales = "free_y") +
  
  # use discrete x with the original numeric labels shown
  scale_x_discrete(labels = c("2","3","4","5","7","10","15")) +
  
  scale_y_continuous(limits = c(0, 1)) +
  
  geom_text(
    data = coxme_results,
    aes(
      x = Inf, y = 0.1,
      label = paste0("p = ", signif(p.value, 2))
    ),
    hjust = 1.1, vjust = 0, size = 3,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Concentration",
    y = paste0("Survival Probability at ", target_time, " time units"),
    title = "Paired Dose–Response Curves with Cox Mixed-Effects (dish random effect)",
    subtitle = "Each facet shows base and its A-variant; p-values from CoxME model"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p1)


# ============================================================
# Plot one KM per concentration using ggsurvplot defaults
# Variants = dashed, base = solid
# ============================================================

library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(stringr)

df <- read_excel("Axenic.xlsx") %>%
  mutate(
    conc = as.numeric(conc),
    dish = as.factor(dish),
    name = str_trim(as.character(name)),
    base = str_replace(name, "A$", ""),
    is_variant = str_detect(name, "A$"),
    cond = ifelse(is_variant, "variant", "base"),
    group = paste(base, cond, sep = "_")
  )

paired_bases <- df %>%
  group_by(base) %>%
  filter(any(is_variant) & any(!is_variant)) %>%
  pull(base) %>%
  unique()

df <- df %>% filter(base %in% paired_bases)

# Precompute linetypes for all groups
all_groups <- unique(df$group)
linetypes <- setNames(
  ifelse(grepl("variant$", all_groups), "dashed", "solid"),
  all_groups
)

for (cval in sort(unique(df$conc))) {
  
  message("Plotting concentration: ", cval)
  
  d_sub <- df %>% filter(conc == cval)
  if (nrow(d_sub) == 0) next
  
  fit_sub <- survfit(Surv(time, death) ~ group, data = d_sub)
  
  p <- ggsurvplot(
    fit_sub,
    data = d_sub,
    conf.int = TRUE,
    legend = "right",
    title = paste("Concentration =", cval),
    xlab = "Time",
    ylab = "Survival probability",
    ggtheme = theme_bw(base_size = 12),
    
    # *** THIS IS THE CRUCIAL FIX ***
    linetype = "group"
  )
  
  p$plot <- p$plot + scale_linetype_manual(values = linetypes)
  
  print(p)
}

