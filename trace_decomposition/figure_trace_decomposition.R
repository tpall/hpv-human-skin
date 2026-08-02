suppressPackageStartupMessages({library(tidyverse); library(jsonlite); library(cowplot); library(here)})
theme_set(theme_cowplot(font_size = 12))
FLOOR <- 13.2
d <- fromJSON(here("trace_decomposition/trace_rows.json")) %>% as_tibble() %>%
  mutate(group = case_when(
    tier %in% c("control_productive","EV") ~ "Produktiivsed\nkontrollid (tüükad/EV)",
    tier %in% c("unselected_skin","neoplasia","other") & rpm < FLOOR ~ "Alla läve\njälg (nahk/vähk)",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(group)) %>%
  mutate(cat = case_when(L1 >= 3 ~ "Produktiivne (≥3 L1-lugemit)",
                         EARLY > 0 ~ "Varajane (E-geenid, L1 puudub)",
                         TRUE ~ "Geenideta (k-meeri müra)"),
         cat = factor(cat, levels = c("Geenideta (k-meeri müra)","Varajane (E-geenid, L1 puudub)","Produktiivne (≥3 L1-lugemit)")),
         group = factor(group, levels = c("Alla läve\njälg (nahk/vähk)","Produktiivsed\nkontrollid (tüükad/EV)")))
comp <- d %>% count(group, cat) %>% group_by(group) %>% mutate(frac = n/sum(n))
write_csv(comp, here("trace_decomposition/trace_composition.csv"))
cols <- c("Geenideta (k-meeri müra)"="grey70","Varajane (E-geenid, L1 puudub)"="#D95F02","Produktiivne (≥3 L1-lugemit)"="#1B9E77")
pA <- ggplot(comp, aes(group, frac, fill = cat)) +
  geom_col(width = .68, colour = "white", linewidth = .3) +
  geom_text(aes(label = ifelse(frac > .04, scales::percent(frac, 1), "")), position = position_stack(vjust=.5), size=3.6) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0,.02))) +
  labs(x = NULL, y = "Osakaal HPV+ teekidest",
       title = "Alla läve jälg = enamjaolt müra + varajane transkriptsioon",
       subtitle = "Produktiivne (L1) programm puudub jäljes, domineerib kontrollides") +
  theme(legend.position = "right")
ggsave(here("trace_decomposition/figure_trace_decomposition.png"), pA, width = 9, height = 5, dpi = 200, bg = "white")
cat("wrote figure + trace_composition.csv\n"); print(comp)
