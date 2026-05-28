library(ggplot2)
library(dplyr)
library(ggnewscale)

# 1. 构造虚拟数据
set.seed(42)
df <- data.frame(
  ID = c(paste0("GO:00", sample(10000:99999, 15)), paste0("hsa0", sample(4000:6000, 5))),
  Ontology = factor(rep(c("Biological Process", "Cellular Component", "Molecular Function", "KEGG Pathway"), c(6, 4, 5, 5)), 
                    levels = c("Biological Process", "Cellular Component", "Molecular Function", "KEGG Pathway")),
  Bg_Genes = sample(50:450, 20, replace = TRUE),
  Select_Genes = sample(10:150, 20, replace = TRUE),
  Pvalue = runif(20, 2, 25),
  Rich_Factor = runif(20, 0.2, 0.9)
)

df <- df %>%
  arrange(Ontology, desc(Bg_Genes)) %>%
  mutate(x = row_number())

# === 关键修改：将 Y 轴生长改为 X 轴 (圆周) 生长 ===
df <- df %>% mutate(
  # 1. 最内圈：Rich Factor (假设这个依然是从内向外生长)
  y1_min = 10, 
  y1_max = 10 + (Rich_Factor * 12),
  x1_min = x - 0.4,
  x1_max = x + 0.4,
  
  # 2. 中间圈 (紫)：Select Genes
  # Y轴(轨道高度)固定。X轴底边紧贴左边缘 (x-0.5)，向右延伸 (最大延伸 0.9，留一点间隙)
  y2_min = 25, 
  y2_max = 33,
  x2_min = x - 0.5,
  x2_max = (x - 0.5) + (Select_Genes / max(Select_Genes)) * 0.9,
  
  # 3. 外圈 (红)：Bg Genes
  # Y轴(轨道高度)固定。X轴底边紧贴左边缘 (x-0.5)，向右延伸
  y3_min = 36, 
  y3_max = 48,
  x3_min = x - 0.5,
  x3_max = (x - 0.5) + (Bg_Genes / max(Bg_Genes)) * 0.9
)

# 2. 绘图
p <- ggplot(df) +
  # === 背景底色与网格 ===
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 10, ymax = 22), fill = "#f5f5f5", color = "white", linewidth = 0.5) +
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 25, ymax = 33), fill = "#ededed", color = "white", linewidth = 0.5) + 
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 36, ymax = 48), fill = "#e5e5e5", color = "white", linewidth = 0.5) + 

  # --- 1. 最内圈：富集因子 ---
  geom_rect(aes(xmin = x1_min, xmax = x1_max, ymin = y1_min, ymax = y1_max, fill = Ontology)) +
  
  # 外围分类底色环
  geom_rect(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = 52, ymax = 58, fill = Ontology), alpha = 0.4, color = "white", linewidth = 0.8) +
  
  scale_fill_manual(values = c("Biological Process" = "#3288bd", 
                               "Cellular Component" = "#99d594", 
                               "Molecular Function" = "#d53e4f", 
                               "KEGG Pathway" = "#fee08b")) +
  new_scale_fill() +
  
  # --- 2. 中间圈 (紫)：沿着扇形宽度向右生长 ---
  # 使用算好的 x2_min 和 x2_max
  geom_rect(aes(xmin = x2_min, xmax = x2_max, ymin = y2_min, ymax = y2_max), fill = "#ba55d3") +
  # 标签移动到柱子的右侧末端
  geom_text(aes(x = x2_max + 0.08, y = (y2_min + y2_max)/2, label = Select_Genes), size = 2.5, fontface = "bold", color = "#8a2be2") +
  
  # --- 3. 外圈 (红)：沿着扇形宽度向右生长 ---
  # 使用算好的 x3_min 和 x3_max
  geom_rect(aes(xmin = x3_min, xmax = x3_max, ymin = y3_min, ymax = y3_max, fill = Pvalue)) +
  scale_fill_gradient(name = "-log10(Pvalue)", low = "#fee0d2", high = "#a50f15") +
  # 标签移动到柱子的右侧末端
  geom_text(aes(x = x3_max + 0.08, y = (y3_min + y3_max)/2, label = Bg_Genes), size = 2.5, color = "black") +
  
  # --- 4. 最外围 ID 标签 ---
  geom_text(aes(x = x, y = 68, label = ID, 
                angle = ifelse(90 - (360 * (x - 0.5) / nrow(df)) < -90, 
                               270 - (360 * (x - 0.5) / nrow(df)), 
                               90 - (360 * (x - 0.5) / nrow(df)))), 
            size = 2.8, fontface = "bold", color = "black") +
  
  coord_polar(theta = "x", start = 0) +
  ylim(0, 75) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# 保存高质量图
ggsave("Enrichment_Ring_Plot_SideBase.pdf", plot = p, width = 10, height = 10)
print("绘图完成！请查看 Enrichment_Ring_Plot_SideBase.pdf")