###FARM-K (Family-Aware Robustness Model for K)

library(ipADMIXTURE)
library(mclust)
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(stringr)
library(clue)

setwd('PATH_TO_DIR')

PED <- 'data'
fam <- read.table(paste0(PED,'.fam'), sep = ' ', header = F)

sample_info <- fam
names(sample_info)[2] <- "SampleID"
names(sample_info)[1] <- "Group"

align_and_average_k <- function(q_list) {
  n_runs <- length(q_list)
  if (n_runs == 1) {
	mat <- as.matrix(q_list[[1]])
	mat <- mat / rowSums(mat)
	colnames(mat) <- paste0("Q", seq_len(ncol(mat)))
	return(mat)
  }
  
  # Нормализация строк на сумму 1
  q_list <- lapply(q_list, function(mat) {
	mat <- as.matrix(mat)
	mat / rowSums(mat)
  })
  
  N <- nrow(q_list[[1]])
  K <- ncol(q_list[[1]])
  
  calc_similarity_matrix <- function(matA, matB) {
	sim_mat <- matrix(0, nrow = K, ncol = K)
	for (i in 1:K) {
	  for (j in 1:K) {
		# Метрика сходства CLUMPP
		sim_mat[i, j] <- 1 - (sum(abs(matA[, i] - matB[, j])) / (2 * N))
	  }
	}
	return(sim_mat)
  }
  
  # Поиск эталонного прогона среди отобранных
  pairwise_scores <- matrix(0, nrow = n_runs, ncol = n_runs)
  for (i in 1:(n_runs - 1)) {
	for (j in (i + 1):n_runs) {
	  sim_mat <- calc_similarity_matrix(q_list[[i]], q_list[[j]])
	  assignment <- clue::solve_LSAP(sim_mat, maximum = TRUE)
	  score <- sum(sim_mat[cbind(1:K, assignment)]) / K
	  pairwise_scores[i, j] <- score
	  pairwise_scores[j, i] <- score
	}
  }
  
  best_run_idx <- which.max(rowMeans(pairwise_scores))
  target_mat <- q_list[[best_run_idx]]
  
  # Выравнивание всех прогонов относительно эталона
  sum_mat <- matrix(0, nrow = N, ncol = K)
  for (r in seq_len(n_runs)) {
	if (r == best_run_idx) {
	  sum_mat <- sum_mat + target_mat
	} else {
	  sim_mat <- calc_similarity_matrix(target_mat, q_list[[r]])
	  perm <- as.vector(clue::solve_LSAP(sim_mat, maximum = TRUE))
	  sum_mat <- sum_mat + q_list[[r]][, perm]
	}
  }
  
  # Консенсусное усреднение
  consensus_q <- sum_mat / n_runs
  consensus_q <- consensus_q / rowSums(consensus_q)
  colnames(consensus_q) <- paste0("Q", seq_len(K))
  rownames(consensus_q) <- rownames(q_list[[1]])
  
  return(consensus_q)
}

# Функция поиска стабильных индивидуумов
find_core_individuals <- function(k_target = 5) {
  cls_mat <- cluster_membership_list[[as.character(k_target)]]
  
  # Проверяем равенство кластеров во всех колонках
  is_stable <- apply(cls_mat, 1, function(row) length(unique(row)) == 1)
  rownames(cls_mat) <- sample_info$SampleID
  df_stability <- data.frame(
	SampleID = rownames(cls_mat),
	Group = sample_info$Group[match(rownames(cls_mat), sample_info$SampleID)],
	AssignedCluster = cls_mat[, 1],
	IsCore = is_stable
  )
  
  # Сводка: процент стабильного ядра по группам
  core_summary <- df_stability %>%
	group_by(Group) %>%
	summarise(
	  Total = n(),
	  Core_N = sum(IsCore),
	  Core_Pct = round(100 * sum(IsCore) / n(), 1)
	) %>%
	arrange(desc(Core_Pct))
  
  return(list(details = df_stability, summary = core_summary))
}

# Функция расчета матрицы Жаккара между двумя векторами кластеризации
calc_jaccard_matrix <- function(cls1, cls2) {
  u1 <- sort(unique(cls1))
  u2 <- sort(unique(cls2))
  
  j_mat <- matrix(0, nrow = length(u1), ncol = length(u2),
				  dimnames = list(paste0("C", u1), paste0("C", u2)))
  
  for (i in seq_along(u1)) {
	set1 <- which(cls1 == u1[i])
	for (j in seq_along(u2)) {
	  set2 <- which(cls2 == u2[j])
	  intersection <- length(intersect(set1, set2))
	  union <- length(union(set1, set2))
	  j_mat[i, j] <- if (union == 0) 0 else intersection / union
	}
  }
  return(j_mat)
}


track_cluster_stability_ths <- function(cls_membership_mat, k_val = 5) {
  base_cls <- cls_membership_mat[, "Ths_0.3"]
  base_clusters <- sort(unique(base_cls))
  ths_cols <- colnames(cls_membership_mat)
  
  res_df <- data.frame()
  
  for (c_id in base_clusters) {
	base_set <- which(base_cls == c_id)
	n_base <- length(base_set)
	
	for (col_name in ths_cols) {
	  curr_cls <- cls_membership_mat[, col_name]
	  j_mat <- calc_jaccard_matrix(base_cls, curr_cls)
	  
	  # Находим наилучшее соответствие в текущем пороге
	  best_match_idx <- which.max(j_mat[paste0("C", c_id), ])
	  max_jaccard <- max(j_mat[paste0("C", c_id), ])
	  
	  ths_num <- as.numeric(gsub("Ths_", "", col_name))
	  
	  res_df <- rbind(res_df, data.frame(
		K = k_val,
		Base_Cluster = paste0("Кластер_", c_id),
		Initial_Size = n_base,
		Threshold = ths_num,
		Jaccard = max_jaccard
	  ))
	}
  }
  return(res_df)
}

track_clusters_across_k <- function(cluster_list, target_ths = "Ths_0.3") {
  k_names <- as.numeric(names(cluster_list))
  k_names <- sort(k_names)
  
  cross_k_results <- data.frame()
  
  for (i in 1:(length(k_names) - 1)) {
	k_curr <- as.character(k_names[i])
	k_next <- as.character(k_names[i + 1])
	
	cls_curr <- cluster_list[[k_curr]][, target_ths]
	cls_next <- cluster_list[[k_next]][, target_ths]
	
	j_mat <- calc_jaccard_matrix(cls_curr, cls_next)
	
	for (row_name in rownames(j_mat)) {
	  c_id <- gsub("C", "", row_name)
	  best_match_col <- colnames(j_mat)[which.max(j_mat[row_name, ])]
	  max_j <- max(j_mat[row_name, ])
	  
	  cross_k_results <- rbind(cross_k_results, data.frame(
		Transition = paste0("K=", k_curr, " -> K=", k_next),
		Source_K = k_curr,
		Target_K = k_next,
		Cluster_ID = paste0("K", k_curr, "_", row_name),
		Size = sum(cls_curr == as.numeric(c_id)),
		Best_Match = paste0("K", k_next, "_", best_match_col),
		Jaccard_Overlap = round(max_j, 3),
		Status = case_when(
		  max_j >= 0.85 ~ "Стабильный",
		  max_j >= 0.50 ~ "Дробящийся/Смешанный",
		  TRUE ~ "Распадающийся/Артефактный"
		)
	  ))
	}
  }
  return(cross_k_results)
}

extract_immutable_clusters <- function(stab_df, threshold_jaccard = 0.85) {
  stab_df %>%
	group_by(K, Base_Cluster) %>%
	summarise(
	  Min_Jaccard = min(Jaccard),
	  Mean_Jaccard = mean(Jaccard),
	  Initial_Size = first(Initial_Size),
	  Is_Robust = all(Jaccard >= threshold_jaccard),
	  .groups = "drop"
	) %>%
	arrange(K, desc(Is_Robust), desc(Initial_Size))
}

# Функция с защитой от пустых совпадений
get_robust_samples_by_k <- function(k_val, cluster_membership_list, robust_summary, sample_info) {
  k_str <- as.character(k_val)
  cls_mat <- cluster_membership_list[[k_str]]
  base_clusters <- cls_mat[, "Ths_0.3"]
  rownames(cls_mat) <- sample_info$SampleID
  
  # Номера стабильных кластеров для этого K
  robust_ids <- robust_summary %>%
	filter(K == k_val & Is_Robust == TRUE) %>%
	mutate(Cluster_Num = as.numeric(gsub("Кластер_", "", Base_Cluster))) %>%
	pull(Cluster_Num)
  
  # Если для данного K нет робастных кластеров, возвращаем пустой data.frame
  if (length(robust_ids) == 0) {
	return(data.frame(
	  SampleID = character(0),
	  K = numeric(0),
	  Cluster = character(0),
	  Group = character(0),
	  stringsAsFactors = FALSE
	))
  }
  
  # Находим индексы образцов
  matched_indices <- which(base_clusters %in% robust_ids)
  
  if (length(matched_indices) == 0) {
	return(data.frame(
	  SampleID = character(0),
	  K = numeric(0),
	  Cluster = character(0),
	  Group = character(0),
	  stringsAsFactors = FALSE
	))
  }
  
  sample_ids <- rownames(cls_mat)[matched_indices]
  
  data.frame(
	SampleID = sample_ids,
	K = k_val,
	Cluster = paste0("Кластер_", base_clusters[matched_indices]),
	Group = sample_info$Group[match(sample_ids, sample_info$SampleID)],
	stringsAsFactors = FALSE
  )
}


#### Подготовка результатов faststructure и pong ####
##### Настоятельно рекомендуется использовать результаты pong
# 1. Загрузка метаданных и лога
summary_lines <- readLines("pong/result_summary.txt") # или result_summary_2.txt
filemap <- read.table("pong_filemap.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(filemap) <- c("RunID", "K", "FilePath")
filemap$Seed <- str_extract(filemap$FilePath, "(?<=faststr_)[0-9]+")

# 2. Парсинг мажорных прогонов и построение выровненного консенсуса
q_consensus_list <- list()
k_indices <- grep("^_{10,}K=", summary_lines)

for (i in seq_along(k_indices)) {
  start_line <- k_indices[i]
  end_line <- if (i < length(k_indices)) k_indices[i + 1] - 1 else length(summary_lines)
  block <- summary_lines[start_line:end_line]
  
  k_val <- as.numeric(str_extract(block[1], "[0-9]+"))
  if (k_val == 1) next
  
  # Определение мажорного режима и его прогонов
  major_line <- grep("^Major mode:", block, value = TRUE)
  major_mode_id <- str_trim(str_remove(major_line, "^Major mode:"))
  
  pattern <- paste0("^\\s*", major_mode_id, "\\s+represents\\s+[0-9]+\\s+runs:\\s+runs\\s+")
  rep_line <- grep(pattern, block, value = TRUE)
  runs_str <- str_remove(rep_line, pattern)
  runs_str <- str_remove(runs_str, "\\..*$")
  major_runs <- str_trim(unlist(str_split(runs_str, ",")))
  
  # Читаем файлы только для валидных мажорных прогонов
  sub_map <- filemap[filemap$RunID %in% major_runs & filemap$K == k_val, ]
  
  q_matrices <- lapply(sub_map$FilePath, function(path) {
	as.matrix(read.table(path, header = FALSE))
  })
  
  cat(sprintf("Выравнивание и усреднение для K = %-2d (%d мажорных прогонов)...\n", 
			  k_val, length(q_matrices)))
  
  # Корректное выравнивание столбцов перед усреднением
  q_consensus_list[[as.character(k_val)]] <- align_and_average_k(q_matrices)
}

# 3. Проверка результата
names(q_consensus_list)
lapply(q_consensus_list, dim)



#### Основной анализ ####
# 1. Задаем сетку параметров
k_values <- 2:14 # Пределы K
ths_values <- seq(0.30, 0.00, by = -0.01) 

# Результирующие структуры
results_grid <- list()
n_clusters_df <- data.frame()
cluster_membership_list <- list()

# 2. Вложенный цикл расчетов
for (k in k_values) {
  # Используем консенсусную Q-матрицу для текущего K
  Q_mat <- as.matrix(q_consensus_list[[as.character(k)]])
  storage.mode(Q_mat) <- "numeric"
  
  cluster_membership_list[[as.character(k)]] <- matrix(
	NA, nrow = nrow(Q_mat), ncol = length(ths_values),
	dimnames = list(rownames(Q_mat), paste0("Ths_", ths_values))
  )
  
  for (t_idx in seq_along(ths_values)) {
	ths <- ths_values[t_idx]
	
	# Запуск иерархической кластеризации
	res <- tryCatch({
	  ipADMIXTURE(
		Qmat = Q_mat, 
		admixRatioThs = ths, 
		method = "average"
	  )
	}, error = function(e) NULL)
	
	if (!is.null(res)) {
	  cls_vec <- res$indexClsVec
	  n_cls <- length(unique(cls_vec))
	  
	  # Сохраняем число кластеров
	  n_clusters_df <- rbind(n_clusters_df, data.frame(
		K = factor(k),
		K_num = k,
		Threshold = ths,
		N_Clusters = n_cls
	  ))
	  
	  # Сохраняем вектор кластеров для оценки стабильности
	  cluster_membership_list[[as.character(k)]][, t_idx] <- cls_vec
	}
  }
}

# 3. Оценка стабильности кластеров
stability_summary <- data.frame()

for (k in k_values) {
  cls_mat <- cluster_membership_list[[as.character(k)]]
  base_clustering <- cls_mat[, "Ths_0.3"]
  
  for (ths_name in colnames(cls_mat)) {
	current_clustering <- cls_mat[, ths_name]
	ari_val <- adjustedRandIndex(base_clustering, current_clustering)
	
	ths_num <- as.numeric(gsub("Ths_", "", ths_name))
	stability_summary <- rbind(stability_summary, data.frame(
	  K = factor(k),
	  K_num = k,
	  Threshold = ths_num,
	  ARI = ari_val
	))
  }
}


# 4. Орисовка графиков
p1 <- ggplot(n_clusters_df, aes(x = Threshold, y = N_Clusters, group = K, color = K)) +
  geom_line(linewidth = .5) +
  geom_point(size = 2) +
  scale_x_reverse(breaks = ths_values) +
  scale_color_viridis_d(option = "plasma", name = "FastSTRUCTURE K") +
  labs(
	title = "Динамика числа кластеров ipADMIXTURE при ужесточении порога примеси",
	subtitle = "Перегибы кривых отражают переход от макропородного разделения к семейному переразбиению",
	x = "Порог доли адмиксии (admixRatioThs)",
	y = "Итоговое число сформированных кластеров"
  ) +
  theme_bw(base_size = 13) +
  theme(
	legend.position = "right",
	panel.grid.minor = element_blank(),
	plot.title = element_text(face = "bold")
  )

print(p1)

p2 <- ggplot(stability_summary, aes(x = factor(Threshold), y = factor(K), fill = ARI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(ARI, 2)), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0, 1), name = "Индекс\nARI") +
  labs(
	title = "Оценка устойчивости предковых кластеров (Adjusted Rand Index)",
	subtitle = "Эталон сравнения: admixRatioThs = 0.30",
	x = "Порог доли примеси (admixRatioThs)",
	y = "Число компонент K (fastSTRUCTURE)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
	panel.grid = element_blank(),
	axis.text = element_text(color = "black", face = "bold"),
	plot.title = element_text(face = "bold")
  )

print(p2)

# Пример анализа выявления ядерных образцов для K = 5
k <- 5
core_k5 <- find_core_individuals(k_target = k)
head(core_k5$summary, 10)

#### Межмодельная стабильность ####

# Пример запуска для K = 5
k <- 5
stab_k5_df <- track_cluster_stability_ths(cluster_membership_list[[as.character(k)]], k_val = k)

# График сохранения кластеров при падении порога
ggplot(stab_k5_df, aes(x = Threshold, y = Jaccard, color = Base_Cluster, group = Base_Cluster)) +
  geom_line(linewidth = .5) +
  geom_point(size = 2.) +
  scale_x_reverse(breaks = seq(0.3, 0.0, -0.05)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  labs(
	title = paste0("Устойчивость состава кластеров при ужесточении порога чистопородности (K=",k,")"),
	subtitle = "Линии выше 0.8 отражают неизменяемые стабильные кластеры",
	color = "Кластер",
	x = "Порог примеси (admixRatioThs)",
	y = "Индекс Жаккара (сходство с базовым составом)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")



#### Cквозноq анализ по всем K
cross_k_stability <- track_clusters_across_k(cluster_membership_list, target_ths = "Ths_0.3")
head(cross_k_stability, 15)
write.table(cross_k_stability, 'межмодельная_стабильность.txt',
			col.names = T, row.names = F, quote = F, sep = '\t')



# 1. Сводный отчет по всем K
all_k_stab_list <- lapply(names(cluster_membership_list), function(k) {
  track_cluster_stability_ths(cluster_membership_list[[k]], k_val = as.numeric(k))
})
all_k_stab_df <- do.call(rbind, all_k_stab_list)

robust_summary <- extract_immutable_clusters(all_k_stab_df, threshold_jaccard = 0.85)
print(robust_summary)

write.table(robust_summary, 'стабильные_кластеры_сводные_результаты.txt',
			col.names = T, row.names = F, quote = F, sep = '\t')


# 2. Сборка общего датафрейма для всех K
robust_samples_list <- lapply(unique(robust_summary$K), function(k) {
  get_robust_samples_by_k(k, cluster_membership_list, robust_summary, sample_info)
})

# Объединяем, отфильтровывая пустые таблицы
robust_samples_all_k <- do.call(rbind, robust_samples_list)

# 3. Проверка результата
head(robust_samples_all_k)
table(robust_samples_all_k$K, robust_samples_all_k$Group)

write.table(robust_samples_all_k, 'список_всех_образцов_стабильных_для_клатера.txt',
			col.names = T, row.names = F, quote = F, sep = '\t')



