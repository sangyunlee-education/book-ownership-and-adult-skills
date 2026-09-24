# 0. 실행 옵션

DATA_PATH   <- "prgkorp2.csv"
USE_VIEWER  <- interactive()
SAVE_OUTPUT <- FALSE
OUTPUT_DIR  <- "output"

# 1. 패키지

packages <- c("tidyverse", "survey", "sensemakr")

to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

# 2. 분석 설정

pv_lit      <- paste0("PVLIT", 1:10)
main_weight <- "SPFWT0"
rep_weights <- paste0("SPFWT", 1:80)
fay_rho     <- 0.3

book_var <- "J2_Q06"
age_var  <- "AGEG5LFS"

book_levels <- c(
  "10권 이하", "11-25권", "26-100권",
  "101-200권", "201-500권", "500권 초과"
)
book_ref <- "10권 이하"

book_col_labels <- c(
  "10권 이하"  = "≤10",
  "11-25권"    = "11-25",
  "26-100권"   = "26-100",
  "101-200권"  = "101-200",
  "201-500권"  = "201-500",
  "500권 초과" = ">500"
)

control_vars <- c(
  "gender_f",
  "age_f",
  "parent_edu_f",
  "mother_work_f",
  "father_work_f",
  "residence_f",
  "family14_f"
)

SENSE_TREAT_TERMS   <- paste0("book_f", setdiff(book_levels, book_ref))
SENSE_BENCHMARK_COV <- "parent_edu_f대졸 이상"
SENSE_KD            <- c(1, 2, 3)
SENSE_ALPHA         <- 0.05

scenario_levels <- paste0(SENSE_KD, "배")

# 3. 보조 함수

as_code <- function(x) suppressWarnings(as.integer(as.character(x)))

recode_char_missing <- function(x) {
  if (is.character(x)) x[x %in% c(".", ".n", ".r", ".d", ".v", "", "NA")] <- NA_character_
  x
}

recode_numeric_special_missing <- function(x) {
  x_num <- as_code(x)
  x_num[x_num %in% c(7, 8, 9, 96, 97, 98, 99, 996, 997, 998, 999)] <- NA_integer_
  x_num
}

make_design <- function(data) {
  survey::svrepdesign(
    weights          = as.formula(paste0("~", main_weight)),
    repweights       = data[, rep_weights],
    data             = data,
    type             = "Fay",
    rho              = fay_rho,
    combined.weights = TRUE
  )
}

make_reg_formula <- function(outcome, controls = control_vars) {
  as.formula(paste(outcome, "~ book_f +", paste(controls, collapse = " + ")))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_dec <- function(x, digits = 3) {
  ifelse(is.na(x), "", sub("^0", "", sprintf(paste0("%.", digits, "f"), x)))
}

p_format <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "<.001",
    TRUE      ~ fmt_dec(p, 3)
  )
}

stars_from_p <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

fmt_p_with_stars <- function(p) paste0(p_format(p), stars_from_p(p))

fmt_ci <- function(low, high) paste0("[", fmt_num(low, 2), ", ", fmt_num(high, 2), "]")

show_table <- function(x, table_name, file_name = NULL) {
  cat("\n\n", table_name, "\n", sep = "")
  print(x, n = Inf, width = Inf)
  if (isTRUE(USE_VIEWER)) View(x, title = table_name)
  if (isTRUE(SAVE_OUTPUT) && !is.null(file_name)) {
    dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
    readr::write_excel_csv(x, file.path(OUTPUT_DIR, file_name))
  }
  invisible(x)
}

# 4. Rubin 결합

pool_svyglm_models <- function(model_list, df_com, outcome_name = "문해력", n_used) {

  coef_names <- Reduce(intersect, lapply(model_list, function(x) names(coef(x))))
  if (length(coef_names) == 0) {
    stop("결합할 공통 회귀계수가 없습니다. 범주 수준 또는 모형식을 확인하세요.")
  }

  coef_mat <- do.call(rbind, lapply(model_list, function(x) coef(x)[coef_names]))
  colnames(coef_mat) <- coef_names

  vcov_list <- lapply(model_list, function(x) {
    as.matrix(vcov(x))[coef_names, coef_names, drop = FALSE]
  })

  k     <- length(model_list)
  q_bar <- colMeans(coef_mat)
  u_bar <- Reduce("+", vcov_list) / k

  if (k > 1) {
    b_mat <- stats::cov(coef_mat)
    if (is.null(dim(b_mat))) {
      b_mat <- matrix(b_mat, nrow = 1, dimnames = list(coef_names, coef_names))
    }
  } else {
    b_mat <- matrix(0, length(q_bar), length(q_bar),
                    dimnames = list(coef_names, coef_names))
  }

  total_var_mat <- u_bar + (1 + 1 / k) * b_mat
  std_error     <- sqrt(diag(total_var_mat))

  b_diag <- diag(b_mat)
  lambda <- pmin(pmax(((1 + 1 / k) * b_diag) / diag(total_var_mat), 0), 0.999999)

  df_old    <- (k - 1) / (lambda^2)
  df_obs    <- ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda)
  df_pooled <- 1 / ((1 / df_old) + (1 / df_obs))
  df_pooled[is.na(df_pooled) | is.infinite(df_pooled)] <- df_com
  df_pooled[b_diag < 1e-12] <- df_com

  statistic <- q_bar / std_error

  tibble(
    outcome        = outcome_name,
    n              = n_used,
    term           = names(q_bar),
    estimate       = as.numeric(q_bar),
    std_error      = as.numeric(std_error),
    conf_low       = as.numeric(q_bar - qt(0.975, df = df_pooled) * std_error),
    conf_high      = as.numeric(q_bar + qt(0.975, df = df_pooled) * std_error),
    statistic_type = "t",
    statistic      = as.numeric(statistic),
    df             = as.numeric(df_pooled),
    p_value        = as.numeric(2 * pt(-abs(statistic), df = df_pooled)),
    stars          = stars_from_p(p_value)
  )
}

pool_scalar_df <- function(df, key_cols, estimate_col, se_col,
                           df_com, outcome_name = "문해력", n_used) {

  stopifnot(length(key_cols) >= 1)

  df %>%
    rename(.est = !!as.name(estimate_col), .se = !!as.name(se_col)) %>%
    filter(!is.na(.est), !is.na(.se)) %>%
    group_by(across(all_of(key_cols))) %>%
    summarise(
      k        = n(),
      estimate = mean(.est),
      u_bar    = mean(.se^2),
      b_var    = if (n() > 1) stats::var(.est) else 0,
      .groups  = "drop"
    ) %>%
    mutate(
      total_var      = u_bar + (1 + 1 / k) * b_var,
      std_error      = sqrt(total_var),
      lambda         = pmin(pmax(((1 + 1 / k) * b_var) / total_var, 0), 0.999999),
      df_old         = (k - 1) / (lambda^2),
      df_obs         = ((df_com + 1) / (df_com + 3)) * df_com * (1 - lambda),
      df             = 1 / ((1 / df_old) + (1 / df_obs)),
      df             = if_else(is.na(df) | is.infinite(df) | b_var < 1e-12, df_com, df),
      conf_low       = estimate - qt(0.975, df = df) * std_error,
      conf_high      = estimate + qt(0.975, df = df) * std_error,
      statistic_type = "t",
      statistic      = estimate / std_error,
      p_value        = 2 * pt(-abs(statistic), df = df),
      stars          = stars_from_p(p_value),
      outcome        = outcome_name,
      n              = n_used
    ) %>%
    select(
      outcome, n, all_of(key_cols),
      estimate, std_error, conf_low, conf_high,
      statistic_type, statistic, df, p_value, stars
    )
}

# 5. 자료 불러오기 및 전처리

read_piaac <- function(path) {
  if (!file.exists(path)) {
    stop("DATA_PATH에 지정한 파일이 없습니다: ", path)
  }

  first_line <- readLines(path, n = 1, warn = FALSE)
  cands      <- c(";", ",", "\t", "|")
  counts     <- vapply(cands, function(d) {
    length(gregexpr(d, first_line, fixed = TRUE)[[1]][
      gregexpr(d, first_line, fixed = TRUE)[[1]] > 0])
  }, integer(1))

  if (max(counts) == 0) {
    stop("구분자를 판별하지 못했습니다. 파일 첫 줄을 확인하세요.")
  }
  delim <- cands[which.max(counts)]

  cat(sprintf("[자료] 구분자 '%s'로 읽습니다.\n",
              if (delim == "\t") "\\t" else delim))

  readr::read_delim(
    path,
    delim          = delim,
    show_col_types = FALSE,
    progress       = FALSE,
    locale         = readr::locale(encoding = "UTF-8")
  )
}

dat_raw <- read_piaac(DATA_PATH)

cat(sprintf("[자료] %s행 %s열\n",
            format(nrow(dat_raw), big.mark = ","),
            format(ncol(dat_raw), big.mark = ",")))

required_raw_vars <- c(
  pv_lit, main_weight, rep_weights, book_var, age_var,
  "GENDER_R", "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

missing_raw_vars <- setdiff(required_raw_vars, names(dat_raw))
if (length(missing_raw_vars) > 0) {
  stop("원자료에서 필요한 변수를 찾지 못했습니다: ",
       paste(missing_raw_vars, collapse = ", "))
}

special_numeric_vars <- c(
  book_var, "PAREDC2", "J2_Q04d", "J2_Q05d",
  "J2_Q07_C", "J2_Q0801", "J2_Q0802"
)

dat <- dat_raw %>%
  mutate(across(where(is.character), recode_char_missing)) %>%
  mutate(across(all_of(special_numeric_vars), recode_numeric_special_missing)) %>%
  mutate(across(all_of(c(pv_lit, main_weight, rep_weights)),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(
    book_f = factor(
      case_when(
        as_code(.data[[book_var]]) == 1 ~ "10권 이하",
        as_code(.data[[book_var]]) == 2 ~ "11-25권",
        as_code(.data[[book_var]]) == 3 ~ "26-100권",
        as_code(.data[[book_var]]) == 4 ~ "101-200권",
        as_code(.data[[book_var]]) == 5 ~ "201-500권",
        as_code(.data[[book_var]]) == 6 ~ "500권 초과",
        TRUE ~ NA_character_
      ),
      levels = book_levels
    ),

    gender_f = factor(
      case_when(
        as_code(GENDER_R) == 1 ~ "남성",
        as_code(GENDER_R) == 2 ~ "여성",
        TRUE ~ NA_character_
      ),
      levels = c("남성", "여성")
    ),

    age_f = factor(
      case_when(
        as_code(.data[[age_var]]) == 1  ~ "20-24세",
        as_code(.data[[age_var]]) == 2  ~ "25-29세",
        as_code(.data[[age_var]]) == 3  ~ "30-34세",
        as_code(.data[[age_var]]) == 4  ~ "35-39세",
        as_code(.data[[age_var]]) == 5  ~ "40-44세",
        as_code(.data[[age_var]]) == 6  ~ "45-49세",
        as_code(.data[[age_var]]) == 7  ~ "50-54세",
        as_code(.data[[age_var]]) == 8  ~ "55-59세",
        as_code(.data[[age_var]]) == 9  ~ "60-64세",
        as_code(.data[[age_var]]) == 10 ~ "65세 이상",
        TRUE ~ NA_character_
      ),
      levels = c("20-24세", "25-29세", "30-34세", "35-39세", "40-44세",
                 "45-49세", "50-54세", "55-59세", "60-64세", "65세 이상")
    ),

    parent_edu_f = factor(
      case_when(
        as_code(PAREDC2) == 1 ~ "중졸 이하",
        as_code(PAREDC2) == 2 ~ "고졸/전문대졸",
        as_code(PAREDC2) == 3 ~ "대졸 이상",
        TRUE ~ NA_character_
      ),
      levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")
    ),

    mother_work_f = factor(
      case_when(
        as_code(J2_Q04d) == 1 ~ "유급직",
        as_code(J2_Q04d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    father_work_f = factor(
      case_when(
        as_code(J2_Q05d) == 1 ~ "유급직",
        as_code(J2_Q05d) == 2 ~ "무직/가사 등",
        TRUE ~ NA_character_
      ),
      levels = c("유급직", "무직/가사 등")
    ),

    residence_f = factor(
      case_when(
        as_code(J2_Q07_C) == 1 ~ "대도시",
        as_code(J2_Q07_C) == 2 ~ "중소도시",
        as_code(J2_Q07_C) == 3 ~ "소도시/읍면",
        as_code(J2_Q07_C) == 4 ~ "농어촌/시골",
        TRUE ~ NA_character_
      ),
      levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")
    ),

    family14_f = factor(
      case_when(
        as_code(J2_Q0801) == 1 & as_code(J2_Q0802) == 1 ~ "양부모 동거",
        as_code(J2_Q0801) == 1 & as_code(J2_Q0802) != 1 ~ "모만 동거",
        as_code(J2_Q0801) != 1 & as_code(J2_Q0802) == 1 ~ "부만 동거",
        as_code(J2_Q0801) != 1 & as_code(J2_Q0802) != 1 ~ "생부모 비동거",
        TRUE ~ NA_character_
      ),
      levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거")
    )
  )

missing_by_var <- dat %>%
  select(all_of(c("book_f", control_vars, pv_lit))) %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "변수", values_to = "결측") %>%
  filter(결측 > 0) %>%
  arrange(desc(결측))

cat("\n[변수별 결측]\n")
print(missing_by_var, n = Inf)

cc_data <- dat %>%
  select(all_of(c(pv_lit, "book_f", control_vars, main_weight, rep_weights))) %>%
  tidyr::drop_na() %>%
  mutate(across(where(is.factor), droplevels)) %>%
  mutate(
    book_f       = relevel(book_f, ref = book_ref),
    parent_edu_f = relevel(parent_edu_f, ref = "중졸 이하"),
    w_           = .data[[main_weight]]
  )

cc_n   <- nrow(cc_data)
cc_des <- make_design(cc_data)
cc_df  <- survey::degf(cc_des)

cat("\n[표본]\n")
cat(sprintf("  원자료           : %s명\n", format(nrow(dat_raw), big.mark = ",")))
cat(sprintf("  결측 제외        : %s명(%.2f%%)\n",
            format(nrow(dat_raw) - cc_n, big.mark = ","),
            100 * (nrow(dat_raw) - cc_n) / nrow(dat_raw)))
cat(sprintf("  최종 분석 표본   : %s명\n", format(cc_n, big.mark = ",")))
cat(sprintf("  복합표본 설계 자유도: %d\n", cc_df))

# 6. 회귀모형 적합 및 결합

main_model_list <- map(pv_lit, function(y) {
  survey::svyglm(make_reg_formula(y), design = cc_des)
})
names(main_model_list) <- pv_lit

main_pooled <- pool_svyglm_models(main_model_list, df_com = cc_df, n_used = cc_n)

main_effect_results <- main_pooled %>%
  filter(term %in% SENSE_TREAT_TERMS) %>%
  mutate(`도서 보유 수` = factor(str_remove(term, "^book_f"),
                            levels = setdiff(book_levels, book_ref))) %>%
  arrange(`도서 보유 수`)

pv_var_vec <- vapply(pv_lit, function(y) {
  as.numeric(survey::svyvar(as.formula(paste0("~", y)), design = cc_des))
}, numeric(1))

sd_lit <- sqrt(mean(pv_var_vec))

cat(sprintf("\n[문해력 가중 표준편차] %.2f점\n", sd_lit))

# 7. 기술통계 및 교차분석

background_specs <- list(
  list(var = "gender_f",      label = "성별",
       levels = c("남성", "여성")),
  list(var = "age_f",         label = "연령",
       levels = c("20-24세", "25-29세", "30-34세", "35-39세", "40-44세",
                  "45-49세", "50-54세", "55-59세", "60-64세", "65세 이상")),
  list(var = "parent_edu_f",  label = "부모 교육수준",
       levels = c("중졸 이하", "고졸/전문대졸", "대졸 이상")),
  list(var = "mother_work_f", label = "모 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "father_work_f", label = "부 경제활동",
       levels = c("유급직", "무직/가사 등")),
  list(var = "residence_f",   label = "거주지역 규모",
       levels = c("대도시", "중소도시", "소도시/읍면", "농어촌/시골")),
  list(var = "family14_f",    label = "가족구조",
       levels = c("양부모 동거", "모만 동거", "부만 동거", "생부모 비동거"))
)

make_row <- function(label, values_chr, p_text = "") {
  values <- as.list(values_chr[book_levels])
  names(values) <- unname(book_col_labels[book_levels])
  bind_cols(tibble(변수 = label), as_tibble(values), tibble(`pᵃ` = p_text))
}

safe_svy_chisq_p <- function(var) {
  fml <- as.formula(paste("~ book_f +", var))
  for (method in c("F", "adjWald", "Chisq")) {
    p <- tryCatch(
      suppressWarnings(as.numeric(survey::svychisq(fml, cc_des, statistic = method)$p.value)),
      error = function(e) NA_real_
    )
    if (length(p) == 1 && is.finite(p)) return(p)
  }
  NA_real_
}

weighted_col_percent <- function(var, category) {
  tmp <- cc_data %>%
    filter(!is.na(book_f), !is.na(.data[[var]])) %>%
    group_by(book_f, category_tmp = as.character(.data[[var]])) %>%
    summarise(weighted_n = sum(.data[[main_weight]], na.rm = TRUE), .groups = "drop") %>%
    group_by(book_f) %>%
    mutate(percent = 100 * weighted_n / sum(weighted_n, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(category_tmp == category)

  out <- setNames(rep(NA_real_, length(book_levels)), book_levels)
  out[as.character(tmp$book_f)] <- tmp$percent
  out
}

make_table1_section <- function(spec) {
  empty <- setNames(rep("", length(book_levels)), book_levels)

  header <- make_row(spec$label, empty, fmt_p_with_stars(safe_svy_chisq_p(spec$var)))
  rows   <- map_dfr(spec$levels, function(cat_label) {
    make_row(cat_label, fmt_num(weighted_col_percent(spec$var, cat_label), 1))
  })

  bind_rows(header, rows)
}

book_distribution <- cc_data %>%
  group_by(book_f) %>%
  summarise(n = n(), weighted_n = sum(.data[[main_weight]]), .groups = "drop") %>%
  mutate(pct = 100 * weighted_n / sum(weighted_n))

n_vec   <- setNames(as.character(book_distribution$n),
                    as.character(book_distribution$book_f))
pct_vec <- setNames(fmt_num(book_distribution$pct, 1),
                    as.character(book_distribution$book_f))

표1_기술통계_교차분석 <- bind_rows(
  make_row("n(비가중)", n_vec),
  make_row("%(가중)",   pct_vec),
  map_dfr(background_specs, make_table1_section)
)

show_table(
  표1_기술통계_교차분석,
  "<표 1> 14세 시점 가정 내 도서 보유 수 범주별 주요 배경변수의 가중 분포",
  "table1.csv"
)

# 8. 회귀분석 결과

표2_회귀분석 <- main_effect_results %>%
  transmute(
    `도서 보유 수`  = as.character(`도서 보유 수`),
    추정치          = fmt_num(estimate, 2),
    표준오차        = fmt_num(std_error, 2),
    t               = fmt_num(statistic, 2),
    p               = fmt_p_with_stars(p_value),
    `95% 신뢰구간`  = fmt_ci(conf_low, conf_high),
    `효과크기(SD)`  = fmt_num(estimate / sd_lit, 2)
  )

show_table(
  표2_회귀분석,
  "<표 2> 14세 시점 가정 내 도서 보유 수가 성인기 문해력에 미치는 효과 추정 결과",
  "table2.csv"
)

# 9. 민감도 분석

lm_model_list <- map(pv_lit, function(y) {
  stats::lm(make_reg_formula(y), data = cc_data, weights = w_)
})
names(lm_model_list) <- pv_lit

r2_y_benchmark <- mean(vapply(lm_model_list, function(m) {
  as.numeric(sensemakr::partial_r2(m, covariates = SENSE_BENCHMARK_COV))
}, numeric(1)))

mm <- model.matrix(make_reg_formula(pv_lit[1]), data = cc_data)
mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

name_map     <- setNames(paste0("v", seq_len(ncol(mm))), colnames(mm))
mm_df        <- as.data.frame(mm)
names(mm_df) <- unname(name_map)
mm_df$w_     <- cc_data$w_

r2_d_benchmark <- vapply(SENSE_TREAT_TERMS, function(tt) {
  y_nm <- name_map[[tt]]
  x_nm <- setdiff(unname(name_map), y_nm)
  m_d  <- stats::lm(as.formula(paste(y_nm, "~", paste(x_nm, collapse = " + "))),
                    data = mm_df, weights = w_)
  as.numeric(sensemakr::partial_r2(m_d, covariates = name_map[[SENSE_BENCHMARK_COV]]))
}, numeric(1))

기준공변량_부분설명력 <- tibble(
  `도서 보유 수`  = str_remove(SENSE_TREAT_TERMS, "^book_f"),
  `R2(원인변수)`  = fmt_dec(r2_d_benchmark, 3),
  `R2(결과변수)`  = fmt_dec(r2_y_benchmark, 3)
)

show_table(
  기준공변량_부분설명력,
  "[보조] 비교기준 공변량('부모 교육수준: 대졸 이상')의 관찰된 부분설명력",
  "benchmark_r2.csv"
)

sense_direction <- setNames(
  ifelse(main_effect_results$estimate >= 0, 1, -1),
  main_effect_results$term
)

make_sense_cache <- function(i, tt) {
  message("민감도 반복추정: ", tt, " / ", pv_lit[i])
  s <- sense_direction[[tt]]

  survey::withReplicates(
    design = cc_des,
    theta = function(w, data) {
      if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
        stop("분석할 수 없는 반복가중치입니다.")
      }
      data$w_ <- as.numeric(w)
      m <- stats::lm(make_reg_formula(pv_lit[i]),
                     data = data, weights = w_)
      tab <- coef(summary(m))
      required <- c(tt, SENSE_BENCHMARK_COV)
      if (!all(required %in% rownames(tab)) || anyNA(coef(m))) {
        stop("반복표본의 모형이 식별되지 않습니다: ", tt)
      }
      beta <- unname(tab[tt, "Estimate"])
      bias_scale <- unname(tab[tt, "Std. Error"] * sqrt(m$df.residual))
      ry <- as.numeric(sensemakr::partial_r2(
        m, covariates = SENSE_BENCHMARK_COV
      ))

      ddata <- mm_df
      ddata$w_ <- as.numeric(w)
      y_nm <- name_map[[tt]]
      x_nm <- setdiff(unname(name_map), y_nm)
      md <- stats::lm(stats::reformulate(x_nm, response = y_nm),
                      data = ddata, weights = w_)
      if (anyNA(coef(md))) stop("원인변수 회귀의 반복표본 특이행렬")
      rd <- as.numeric(sensemakr::partial_r2(
        md, covariates = name_map[[SENSE_BENCHMARK_COV]]
      ))

      adj <- vapply(SENSE_KD, function(k) {
        b <- sensemakr::ovb_partial_r2_bound(
          r2dxj.x = rd, r2yxj.dx = ry, kd = k, ky = k
        )
        r_d <- b$r2dz.x
        r_y <- b$r2yz.dx
        if (length(r_d) != 1L || length(r_y) != 1L ||
            any(!is.finite(c(r_d, r_y))) ||
            r_d < 0 || r_d >= 1 || r_y < 0 || r_y > 1) {
          stop("비교기준 시나리오가 허용 범위를 벗어났습니다: ", tt)
        }
        beta - s * bias_scale * sqrt(r_y * r_d / (1 - r_d))
      }, numeric(1))

      out <- c(beta = beta, bias_scale = bias_scale,
               setNames(adj, paste0("adj_", seq_along(SENSE_KD))))
      if (any(!is.finite(out))) stop("민감도 계산에 비유한 값이 있습니다.")
      out
    },
    return.replicates = TRUE
  )
}

sense_cache <- setNames(lapply(SENSE_TREAT_TERMS, function(tt) {
  lapply(seq_along(pv_lit), function(i) make_sense_cache(i, tt))
}), SENSE_TREAT_TERMS)

pool_sense <- function(tt, f = 0, j = NULL) {
  s <- sense_direction[[tt]]
  pv_rows <- lapply(seq_along(pv_lit), function(i) {
    obj <- sense_cache[[tt]][[i]]
    z <- coef(obj)
    reps <- as.matrix(obj$replicates)
    if (ncol(reps) != length(z)) stop("반복추정치 열 수 불일치")
    colnames(reps) <- names(z)
    if (is.null(j)) {
      q <- unname(z["beta"] - s * f * z["bias_scale"])
      qr <- reps[, "beta"] - s * f * reps[, "bias_scale"]
    } else {
      nm <- paste0("adj_", j)
      q <- unname(z[nm])
      qr <- reps[, nm]
    }
    if (!is.finite(q) || any(!is.finite(qr))) {
      stop("반복추정치가 누락되었습니다: ", tt)
    }
    u <- as.numeric(survey::svrVar(
      thetas = qr, scale = cc_des$scale, rscales = cc_des$rscales,
      mse = cc_des$mse, coef = q
    ))
    if (!is.finite(u) || u < 0) stop("유효하지 않은 반복분산")
    tibble(treatment = tt, pv = pv_lit[i], est = q, se = sqrt(u))
  })
  ans <- pool_scalar_df(
    df = bind_rows(pv_rows), key_cols = "treatment",
    estimate_col = "est", se_col = "se",
    df_com = cc_df, n_used = cc_n
  )
  crit <- qt(1 - SENSE_ALPHA / 2, df = ans$df)
  ans$conf_low <- ans$estimate - crit * ans$std_error
  ans$conf_high <- ans$estimate + crit * ans$std_error
  ans
}

for (tt in SENSE_TREAT_TERMS) {
  z <- pool_sense(tt, f = 0)
  ref <- main_effect_results[main_effect_results$term == tt, ]
  if (nrow(ref) != 1L ||
      !isTRUE(all.equal(z$estimate, ref$estimate, tolerance = 1e-6)) ||
      !isTRUE(all.equal(z$std_error, ref$std_error, tolerance = 1e-6))) {
    stop("교란 0 재현 점검 실패: ", tt,
         ". 자료 순서와 설계 설정을 확인하세요.")
  }
}
message("교란 0에서 원래 회귀분석 재현: 통과")

sensitivity_bounds_pooled <- bind_rows(lapply(SENSE_TREAT_TERMS, function(tt) {
  bind_rows(lapply(seq_along(SENSE_KD), function(j) {
    z <- pool_sense(tt, j = j)
    z$scenario <- scenario_levels[j]
    z
  }))
}))

find_sense_threshold <- function(tt) {
  s <- sense_direction[[tt]]
  margin <- function(f) {
    z <- pool_sense(tt, f = f)
    s * z$estimate - qt(1 - SENSE_ALPHA / 2, z$df) * z$std_error
  }
  if (margin(0) <= 0) return(0)
  z0 <- pool_sense(tt, f = 0)
  mean_scale <- mean(vapply(sense_cache[[tt]], function(x) {
    unname(coef(x)["bias_scale"])
  }, numeric(1)))
  if (!is.finite(mean_scale) || mean_scale <= 0) stop("편향 척도 오류")
  f_zero <- abs(z0$estimate) / mean_scale
  grid <- seq(0, f_zero * (1 + 1e-8), length.out = 257)
  vals <- vapply(grid, margin, numeric(1))
  hit <- which(vals <= 0)[1]
  if (is.na(hit)) stop("신뢰구간 임계값을 찾지 못했습니다: ", tt)
  uniroot(margin, interval = grid[c(hit - 1L, hit)], tol = 1e-8)$root
}

sensitivity_thresholds <- bind_rows(lapply(SENSE_TREAT_TERMS, function(tt) {
  f_star <- find_sense_threshold(tt)
  tibble(
    treatment = tt,
    r_equal_ci = if (f_star == 0) 0 else
      2 / (1 + sqrt(1 + 4 / f_star^2)),
    r_d_ci_y1 = f_star^2 / (1 + f_star^2)
  )
}))

표3_벤치마크 <- sensitivity_bounds_pooled %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    조정치 = paste0(fmt_num(estimate, 2), stars_from_p(p_value),
                   " (", fmt_num(std_error, 2), ")")
  ) %>%
  select(treatment, scenario, 조정치) %>%
  pivot_wider(names_from = scenario, values_from = 조정치,
              names_glue = "{scenario} 조정치(SE)")

표3_민감도분석 <- sensitivity_thresholds %>%
  left_join(표3_벤치마크, by = "treatment") %>%
  transmute(
    `도서 보유 수` = str_remove(treatment, "^book_f"),
    `동일 부분R2 CI 임계값` = fmt_dec(r_equal_ci, 3),
    `rY=1일 때 rD CI 임계값` = fmt_dec(r_d_ci_y1, 3),
    across(ends_with("조정치(SE)"))
  )

show_table(표3_민감도분석,
           "<표 3> 민감도 조정 통계량의 반복추정 결과", "table3.csv")

# 10. 실행 환경

cat("\n[재현 정보]\n")
cat("  ", R.version.string, "\n", sep = "")
for (pkg in packages) {
  cat(sprintf("   %-12s %s\n", pkg, as.character(packageVersion(pkg))))
}
