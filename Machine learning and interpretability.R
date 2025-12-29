
# 安装R包
install.packages("pbapply")
install.packages("openxlsx")
install.packages("gbm")
install.packages("caret")
install.packages("rms")
install.packages("rmda")
install.packages("dcurves")
install.packages("ResourceSelection")
install.packages("survey")
install.packages("plotroc")
install.packages("shapper")
install.packages("iml")
install.packages("e1071")
install.packages("ROCR")
install.packages("corrplot")
install.packages("lattice")
install.packages("Formula")
install.packages("SparseM")
install.packages("survival")
install.packages("riskRegression")
install.packages("pheatmap")
install.packages("fastshap")
install.packages("table1")
install.packages("tableone")
install.packages("adabag")
install.packages("RColorBrewer")
install.packages("VIM")
install.packages("mice")
install.packages("autoReg")
install.packages("cvms")
install.packages("tibble")
install.packages("plotROC")
install.packages("pROC")
install.packages("data.table")
install.packages("circlize")
install.packages("ROSE")
install.packages("DMwR")
install.packages("scales")
install.packages("kernelshap")
install.packages("xts")
install.packages("quantmod")
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("rpart")
install.packages("rpart.plot")
install.packages("ggplot2")
install.packages("randomForest")
install.packages("xgboost")
install.packages("lightgbm")
install.packages("kknn")
install.packages("neuralnet")
install.packages("NeuralNetTools")
install.packages("gridExtra")
install.packages("partykit")
install.packages("missForest")
install.packages("regplot")


# 加载R包
library(pbapply)
library(rlang)
library(reshape2)
library(openxlsx)
library(DALEX)
library(readr)
library(gbm)
library(dplyr)
library(caret)
library(ggplot2)
library(pROC)
library(rms)
library(rmda)
library(dcurves)
library(Hmisc)
library(ResourceSelection)
library(survey)
library(foreign)
library(plotROC)
library(shapper)
library(iml)
library(e1071)
library(ROCR)
library(corrplot)
library(lattice)
library(Formula)
library(SparseM)
library(survival)
library(riskRegression)
library(pheatmap)
library(fastshap)
library(ingredients)
library(mlr3)
library(table1)
library(tableone)
library(adabag)
library(RColorBrewer)
library(VIM)
library(mice)
library(autoReg)
library(cvms)
library(tibble)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(ROSE)
library(DMwR)
library(scales)
library(kernelshap)
library(shapviz)
library(rpart)       
library(rpart.plot) 
library(randomForest)  
library(xgboost)         
library(lightgbm)      
library(kknn)        
library(neuralnet)    
library(NeuralNetTools) 
library(gridExtra) 
library(partykit)
library(missForest)
library(regplot)
library(glmnet)


# 设置工作空间 (数据存储和结果输出的文件夹)
setwd("D://Machine learning//Machine learning + interpretability")
# 查看工作空间
getwd()   

# 读取数据
log <-read.csv(file="diabetes.csv",header = T,encoding = "GBK")

# 查看导入文件的情况
head(log)  # 查看前六行数据
str(log)   # 查看数据结构

####################1.数据预处理################################################

##1.1 数据插补(随机森林插补)#####

# 检查缺失值
missing_data <- sapply(log, function(x) sum(is.na(x)))  # 计算数据框每一列的缺失值数量
print(missing_data)

# 插补数据
set.seed(123)    # 设置随机种子
missForest(log)  # 填补数据
log_imputed <- missForest(log)$ximp   # 提取插补后的数据

# 查看插补后的变量
log_imputed$insulin 

# sp为整数,将插补后的数据四舍五入为整数
log_imputed$insulin <- round(log_imputed$insulin)
log_imputed$insulin 

# 查看插补后数据的缺失情况
print(sapply(log_imputed, function(x) sum(is.na(x))))

#导出随机森林插补后的数据
write.csv(log_imputed,file = "logistic_imputed.csv",row.names = FALSE)

# 读取插补后的完整数据 
log <-read.csv(file="logistic_imputed.csv",header = T,encoding = "GBK")

##1.2 数据标准化#####
#对自变量中的定量变量进行标准化处理,适用于基于距离度量的knn、svm等模型
names(log)   # 返回数据框中所有列的名称
set.seed(123)  
# 对log数据框中的定量自变量进行标准化
log$age_scaled <- scale(log$age)[, 1]
log$glucose_scaled <- scale(log$glucose)[, 1]
log$pressure_scaled <- scale(log$pressure)[, 1]
log$triceps_scaled <- scale(log$triceps)[, 1]
log$insulin_scaled <- scale(log$insulin)[, 1]
log$bmi_scaled <- scale(log$bmi)[, 1]
log$pedigree_scaled <- scale(log$pedigree)[, 1]

# 查看结果
head(log)

##1.3 数据集拆分#####
set.seed(12) # 设置随机种子
randnum <- createDataPartition(y=log$diabetes, # 指定目标变量,根据该列分层采样
                               p=0.70,  # 划分比例,70%的数据将用于训练集,30%测试集
                               list = FALSE   # 决定函数返回的结果格式,是一个向量
                               )  

#根据随机数字,产生训练集(tlog)
tlog<-log[randnum,]         # 从log数据框中选择索引randnum指定的行

#根据随机数字,产生验证集(valdata)
valdata<-log[-randnum,]     #  提取log中不在randnum索引中的行,剩余的30%数据

# 查看结局变量的分布特征
table(tlog$diabetes )              # 频数分布
prop.table(table(tlog$diabetes ))  # 相对频率

table(valdata$diabetes )         
prop.table(table(valdata$diabetes ))


##1.4 训练集与验证集均衡性比较#####
# 两个数据集都生成一个group变量
tlog$group <- "训练集"
valdata$group <- "验证集"

# 合并数据集,训练集
total <- rbind(tlog, valdata)
names(total)

# 创建训练集与验证集描述性统计表
baseline <- CreateTableOne(vars = c("diabetes","gender","exercise","race","his",
                                    "hyperlip","pregnant","age", "glucose",
                                    "pressure","triceps","insulin","bmi",
                                    "pedigree"), #描述性统计的变量
                         strata = "group",   ## 按照group进行分组
                         data = total,
                         factorVars = c("diabetes","gender","exercise","race",
                                        "his","hyperlip"))  # 作为因子变量处理的变量

# 输出表格结果
baseline_table  <- print(baseline) 
# 导出为 CSV 文件
write.csv(baseline_table,file = '训练集与验证集均衡性比较结果.csv') 


######################2.特征选择##############################################

#基于训练集进行特征选择 (特征选择方法很多,这里给出常用的两种)

##2.1 先单后多#####
univariate <- CreateTableOne(vars = c("gender","exercise","race","his",
                                      "hyperlip","pregnant","age", "glucose",
                                      "pressure","triceps","insulin","bmi",
                                      "pedigree"), #描述性统计的变量
                           strata = "diabetes",   # 按照diabetes进行分组
                           data = tlog,
                           factorVars = c("gender","exercise","race","his",
                                          "hyperlip"))  # 作为因子变量处理的变量

# 输出表格结果
univariate_table  <- print(univariate) 
# 导出为 CSV 文件
write.csv(univariate_table,file = '单因素分析结果.csv') 

# 多因素logistic回归

# 数据处理：分类自变量因子化
tlog$gender1 <- factor(tlog$gender,levels = c(0,1),labels = c('Female','Male')) 
tlog$exercise1 <- factor(tlog$exercise,levels = c(0,1),labels = c('No','Yes')) 
tlog$race1<-factor(tlog$race,levels = c(1,2,3),labels=c('White','Black','Other')) 
tlog$his1 <- factor(tlog$his,levels = c(0,1),labels = c('No','Yes')) 
tlog$hyperlip1 <- factor(tlog$hyperlip,levels = c(0,1),labels = c('No','Yes')) 

# 构建多因素logistic回归模型
multivariate <- glm(diabetes ~ gender1 +exercise1 + race1 + his1 + hyperlip1 +
                     pregnant + age +glucose + pressure + triceps +insulin +
                     bmi + pedigree,   # y ~ x
                   data=tlog,          # 拟合模型的数据集
                   family=binomial     # 使用二项分布，适用于二分类因变量
                   ) 
# 查看回归模型摘要
summary(multivariate)

# 提取p值
p_values <- summary(multivariate)$coefficients[, 4]

# 过滤出 p 值小于 0.05 的变量
significant_vars <- names(p_values)[p_values < 0.05]

# 显示显著变量
print(significant_vars)


## autoReg 函数生成单因素和多因素Logistic结果
LR_table <- autoReg(multivariate,   # 已建立的多因素 logistic 回归模型
                    uni=TRUE,       # 进行单因素回归分析
                    milti=TRUE,     # 进行多因素回归分析
                    threshold=0.05  # 显著性水平阈值0.05
                    )
LR_table
write.csv(LR_table,"先单后多Logistic回归.csv",row.names = F) #保存单、多因素logistic回归结果


##2.2 LASSO#####

# 提取自变量（第2列~第14列）和因变量（第1列变量）
x <- as.matrix(tlog[, 2:14])  # 自变量（矩阵形式）
y <- tlog$diabetes            # 因变量

# 执行 Lasso 回归
lasso_model <- glmnet(x, y, 
                      family = "binomial",  # 表示二分类回归
                      alpha = 1             # alpha=1表示Lasso
                      ) 

# 绘制 Lasso 回归的系数路径图
plot(lasso_model,      #  Lasso 模型
     xvar = "lambda",  # 以lambda(正则化参数)为 x 轴
     label = T         # 在图中标注系数的名称
     )

# 查看 Lasso 模型的系数
print(lasso_model)

# 使用交叉验证来选择最佳的 λ
set.seed(123)  # 设置随机种子
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
cv_lasso 
plot(cv_lasso) # 绘制交叉验证误差图

# 获取最佳 λ 值
best_lambda <- cv_lasso$lambda.min
best_lambda
# 提取系数（包括零系数）
coef_lasso <- coef(cv_lasso, s = best_lambda )   # 可以修改：s = 1se
coef_lasso
# 将系数转换为矩阵
coef_lasso_matrix <- as.matrix(coef_lasso)

# 选择非零系数对应的特征
selected_features <- rownames(coef_lasso_matrix)[          # 获取矩阵的行名（即变量的名称）
                             coef_lasso_matrix[, 1] != 0   # 筛选出非零系数对应的变量
                             ]  

# 显示选择的特征
print(selected_features)
cat("选择的特征为：",selected_features,"\n")


####
# 依据特征选择的变量,构建机器学习模型
# 本次课程依据先单后多结果,提取多因素Logistic有意义的变量(可换成lasso筛选的变量)
# 变量："exercise", "hyperlip","pregnant","age", "glucose", "bmi", "pedigree" 

# 定义筛选特征集合
selected_vars <- c("exercise", "hyperlip", "pregnant", 
                   "age", "glucose", "bmi", "pedigree")  # 原始特征
selected_vars_scaled <- c("exercise", "hyperlip", "pregnant", "age_scaled", 
                "glucose_scaled", "bmi_scaled", "pedigree_scaled")  #标准化特征


# 将结局变量因子化
tlog$diabetes <- factor(tlog$diabetes,levels = c(0,1),labels = c('No','Yes'))
valdata$diabetes <- factor(valdata$diabetes,levels = c(0,1),labels = c('No','Yes'))
### 注意: 构建模型时,自变量中的分类自变量没有因子化,原因是后续shap法无法识别因子变量)


####################3.二分类机器学习模型建模####################################

# 基于训练集构建模型

###### 3.1 Logistic模型 #########

# 拟合模型
lr_model<- glm(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree,
               data = tlog,
               family ="binomial"  #使用二项分布，适用于二分类因变量
               ) 

# 显示模型信息 
print(lr_model)    

## 绘制列线图
regplot(lr_model, 
        title = "Nomogram", 
        points = TRUE,                    # 显示每个变量的点数贡献
        axis.text.size = 12,              # 调整刻度字体大小
        title.text.size = 14)             # 调整标题字体大小

################3.2 决策树:分类回归树########################################

# 构建基础CART模型:利用默认参数建模
tree_model1 <- rpart(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree, 
                     data = tlog, 
                     method = "class")   # 分类问题,使用分类树算法来构建决策树
tree_model1$cptable   # 返回模型剪枝的复杂度表

# 设置控制参数  
control <- trainControl(method="cv", number=10)  
# 定义参数网格
param_grid <- expand.grid(cp = seq(0.001, 0.3, by = 0.002))  #cp:CART模型的复杂度参数,控制模型的剪枝过程

# 使用train函数进行交叉验证和模型调优
set.seed(111)
fit_cv_rpart <- train(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree, 
                      data = tlog, 
                      method = "rpart", # 使用 rpart算法来训练模型
                      trControl = control, 
                      tuneGrid = param_grid)
fit_cv_rpart $bestTune  # 查看最优参数

# 使用最佳参数构建决策树模型  
tree_model <- rpart(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree, 
                    data = tlog, 
                    method = "class", 
                    cp=fit_cv_rpart$bestTune) 

# 查看控制参数
print(tree_model$control)

# 画决策树图
plot(as.party(tree_model))

# 显示模型信息 
print(tree_model)


###################3.3 随机森林(RF)模型##########################

# 构建默认参数,构建基础随机森林模型
rf_model0 <- randomForest(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree, 
                          data = tlog,  
                          importance=TRUE) # 利用默认参数建模
print(rf_model0)

## 最佳模型参数：超参数调节、网格搜索、交叉验证
# 定义训练控制参数
set.seed(123)  
ctrl <- trainControl(method = "cv", 
                     number = 10,    # 10折交叉验证
                     search = "grid")  # 网格搜索

# 定义超参数mtry搜索范围(mtry表示每棵树随机选择的特征数)
tuneGrid <- expand.grid(mtry = c(1:sqrt(7)))   # 从1到数据集中特征数的平方根的整数值        
                                   
# 超参数调优
rf_model1 <- train(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree, 
                  data = tlog, 
                  method = "rf",      # 指定使用 rf
                  trControl = ctrl,   # 指定训练控制参数 ctrl
                  tuneGrid = tuneGrid # 指定要调优的超参数网格
                  )
# 输出最佳模型参数
print(rf_model1)
rf_model1 $ bestTune

# 设置树的数量ntree范围
ntree_values <- seq(50, 1000, by = 50)  
# 创建一个向量,存储每个 ntree 值对应的 OOB 错误率
oob_error_rates <- numeric(length(ntree_values))

# 训练多个模型并记录OOB误差率
for (i in 1:length(ntree_values)) {   # for循环：使用不同的ntree训练多个模型
  rf_model2 <- randomForest(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree,
                           data = tlog, 
                           mtry = rf_model1$bestTune$mtry, 
                           ntree = ntree_values[i], 
                           importance = TRUE, 
                           oob.prox = TRUE)
  # 找出OOB误差率最低的树的数量
  oob_error_rates[i] <- rf_model2$err.rate[ntree_values[i]]
}

# 找出OOB误差率最低的树的数量
best_ntree <- ntree_values[which.min(oob_error_rates)]
print(paste("最佳树的数量：", best_ntree))

# 使用最佳参数重构建模型
rf_model <- randomForest(diabetes ~ exercise+hyperlip+pregnant+age+glucose+bmi+pedigree,
                         data = tlog, 
                         ntree = best_ntree, 
                         mtry = rf_model1$bestTune$mtry,
                         importance = TRUE)

# 显示模型信息 
print(rf_model)


########3.4 Xgboost模型################################

## 因变量需为数值型变量
tlog$diabetes <- as.numeric(tlog$diabetes) - 1
valdata$diabetes <- as.numeric(valdata$diabetes) - 1

# 设置XGBoost的训练和验证数据集 
train_matrix <- xgb.DMatrix(data = as.matrix(tlog[, selected_vars]), 
                             label = tlog$diabetes)  # label指定模型结局变量

val_matrix <- xgb.DMatrix(data = as.matrix(valdata[, selected_vars]), 
                           label = valdata$diabetes)  

# 基于默认参数,构建基础 Xgboost模型
xgb_model0 <- xgb.train(data = train_matrix, nrounds=100) # 100次迭代


# 超参数调优  
param_grid <- expand.grid(
  objective = "binary:logistic",  # 二分类任务,预测类别概率
  max_depth = c(2, 3, 4, 5),      # 树的最大深度,控制模型的复杂度
  eta = c(0.01, 0.1, 0.2),        # 学习率,决定每一轮迭代中模型更新的步伐大小
  nrounds = c(50, 100, 150)       # 训练轮数,每轮会调整树的参数
)  

# 初始化最佳 AUC 和参数
best_auc <- 0  
best_params <- list()     # 存储最优超参数组合

# 超参数调优的循环
for (i in 1:nrow(param_grid)) {  
  param <- list(  
    objective = "binary:logistic",  
    eval_metric = "auc",  
    max_depth = param_grid$max_depth[i],  
    eta = param_grid$eta[i]  
  )  
  
  xgb_model_0 <- xgb.train(params = param, data = train_matrix, 
                           nrounds = param_grid $ nrounds[i])  
  
  # 评估每个模型的 AUC  
  pred_probs <- predict(xgb_model_0, train_matrix)  
  roc_curve <- roc(tlog$diabetes, pred_probs)  
  auc_value <- roc_curve$auc  
  
  #  选择最优模型
  if (auc_value > best_auc) {  
    best_auc <- auc_value  
    best_params <-  c(param, nrounds = param_grid$nrounds[i])  
  }  
} 

# 输出最佳参数和AUC  
print(best_params)  # 输出最优的超参数组合 
cat("最佳AUC: ", best_auc, "\n")  

# 使用最佳超参数训练模型  
xgb_model <- xgb.train(params = best_params, data = train_matrix, 
                       nrounds = best_params$nrounds)  

#显示模型信息 
print(xgb_model)


#######################3.5 LightGBM模型#############################

# 设置LightGBM 的训练和验证数据集 
lgbtlog <- lgb.Dataset(as.matrix(tlog[,selected_vars]),
                       label = tlog$diabetes) # 创建LightGBM所需的训练数据格式
lgbvaldata <- lgb.Dataset.create.valid(lgbtlog, 
                           as.matrix(valdata[,selected_vars]), 
                           label = tlog$diabetes) # 创建验证集,用于模型的验证过程

# 基于默认参数,构建基础 LightGBM 模型
lightgbm_model0 <- lgb.train(data = lgbtlog)  


# 设置超参数搜索网格  
param_grid <- expand.grid(  
  num_leaves = c(15, 31),    # 树的叶子数
  max_depth = c(-1, 1, 3),   # 树的最大深度
  learning_rate = c( 0.1, 0.2),  # 学习率,控制每一轮迭代时模型更新的步伐
  n_estimators = c(50),          # 训练轮数（树的数量）
  min_data_in_leaf = c(30),      # 每棵树叶子节点最少样本数
  lambda_l1 = c(0, 1),           # L1 正则化参数
  lambda_l2 = c(0, 1)            # L2 正则化参数
)  

# 准备一个空的数据框来保存结果  
results <- data.frame()  

# 超参数调优与交叉验证 
for (i in 1:nrow(param_grid)) {  
  
  # 获取当前的参数组合  
  params <- list(  
    objective = "binary",  
    metric = "auc",  
    learning_rate = param_grid$learning_rate[i],  
    num_leaves = param_grid$num_leaves[i],  
    max_depth = param_grid$max_depth[i],  
    n_estimators = param_grid$n_estimators[i],  
    min_data_in_leaf = param_grid$min_data_in_leaf[i]  
  )  
  
  # 进行交叉验证  
  cv_results <- lgb.cv(  
    params = params,  
    data = lgbtlog,  
    nrounds = 10,  
    nfold = 5,  
    early_stopping_rounds = 10,  
    verbose = -1    # -1 表示不输出训练过程的详细信息
  )  
  
  # 保存当前的参数和其对应的auc  
  results <- rbind(results,data.frame(param_grid[i, ], 
                             auc=max(cv_results$record_evals$valid[['auc']]$data)))  
}  

# 找到最优参数  
best_params <- results[which.max(results$auc), ]  
print(best_params)  

# 用最佳参数构建模型  
best_params_list <- list(  #超参数配置列表
  objective = "binary",  
  metric = "auc",  
  learning_rate = best_params$learning_rate,  
  num_leaves = best_params$num_leaves,  
  max_depth = best_params$max_depth,  
  n_estimators = best_params$n_estimators,  
  min_data_in_leaf = best_params$min_data_in_leaf  
)  

lightgbm_model <- lgb.train(  
  params = best_params_list,  
  data = lgbtlog,  
  nrounds = best_params$n_estimators  
)  

#显示模型信息 
print(lightgbm_model)


########################3.6 knn 模型###########################

# 将结局变量因子化
tlog$diabetes <- factor(tlog$diabetes,levels = c(0,1),labels = c('No','Yes'))

# 基于默认参数,构建基础 knn 模型
knn_model0 <- train(diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                     glucose_scaled + bmi_scaled + pedigree_scaled,  
                   data = tlog, 
                   method = "kknn"    # 指定使用加权K近邻算法
                   )


# 设置交叉验证控制
train_control <- trainControl(method = "cv", number = 5)

# 设置超参数网格，核函数和 k 值
tune_grid<-expand.grid(kmax = seq(1, 20, by = 2), # 调整 k 值
                  distance = 2,             # Minkowski距离，2表示欧几里得距离
                  kernel=c("rectangular","triangular","gaussian")) # 核函数,计算邻居的权重

# 训练 KNN 模型并调优
set.seed(123)
kknn_model <- train(diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                      glucose_scaled + bmi_scaled + pedigree_scaled,  
                    data = tlog, 
                    method = "kknn", 
                    trControl = train_control, 
                    tuneGrid = tune_grid)

# 查看调参结果
print(kknn_model)

# 提取最佳参数组合
best_params <- kknn_model$bestTune
print(best_params)

# 绘制可视化调参结果
ggplot(kknn_model) +
  theme_minimal() +
  ggtitle("KNN 超参数调整结果")

# 使用最佳参数构建最终模型
knn_model <- train(diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                     glucose_scaled + bmi_scaled + pedigree_scaled, 
                   data = tlog,  
                   method = "kknn", 
                   trControl = train_control, 
                   tuneGrid = expand.grid(kmax = best_params$kmax, 
                                          distance = best_params$distance, 
                                          kernel = best_params$kernel))  

# 查看最终模型
print(knn_model)


#########################3.7 支持向量机(SVM)#######################################################
##基于标准化后的数据建模

# 参数调整
set.seed(11)
tune_result <- tune.svm(diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                          glucose_scaled + bmi_scaled + pedigree_scaled, 
                        data = tlog,   
                        kernel = "radial",   # 径向核函数RBF
                        cost = 10^(-1:3),  # cost:惩罚参数,用于控制分类错误的惩罚程度
                        gamma = 10^(-3:1), # gamma:核函数的参数,定义单个训练样本影响的范围
                        tunecontrol=tune.control(sampling = "cross",cross = 5), #交叉验证
                        probability = TRUE) 

# 查看最佳参数  
best_model <- tune_result$best.model  
print(tune_result)  

# 使用最佳参数拟合SVM模型  
svm_model <- svm(diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                   glucose_scaled + bmi_scaled + pedigree_scaled, 
                 data = tlog,   
                 kernel = "radial",   
                 cost = best_model$cost,   
                 gamma = best_model$gamma,   
                 probability = TRUE)  # 启用概率 
print(svm_model)


#########################3.8 神经网络(nnet)#######################################################

# 构建神经网络模型的函数  
build_nn_model <- function(hidden_layers) {  
  formula <- as.formula("diabetes ~ exercise + hyperlip + pregnant + age_scaled + 
                          glucose_scaled + bmi_scaled + pedigree_scaled")  
  model <- neuralnet(formula, data = tlog, hidden = hidden_layers, 
                     linear.output = FALSE)  # 模型输出为分类概率(非线性激活)
  return(model)  
}  

# 初始化变量  
best_model_nnet <- NULL  
best_auc <- 0  
best_hidden_layers <- NULL  # 保存最佳隐藏层组合

# 设置隐藏层组合  
hidden_layer_combinations <- list(c(2),c(3),c(4), c(2, 1))  
# 网格搜索
for (hidden in hidden_layer_combinations) {  # 遍历每种隐藏层结构
  set.seed(123)  # 设置随机种子
  nn_model <- build_nn_model(hidden)  
  
  # 进行预测
  predictions_prob <- predict(nn_model, tlog)[,2]   # 获取概率  
  predictions <- ifelse(predictions_prob > 0.5, "Yes", "No")  # 将概率转为分类  
  
  # 计算AUC  
  roc_obj <- roc(tlog$diabetes, predictions_prob)  
  auc_value <- roc_obj$auc  
  
  # 更新最佳模型  
  if (auc_value > best_auc) {  
    best_auc <- auc_value  
    best_model_nnet <- nn_model
    best_hidden_layers <- hidden  # 保存最佳隐藏层组合
  }  
}  

# 输出最佳模型和AUC值  
cat("Best AUC:", best_auc, "\n")
cat("Best Hidden Layer Configuration:", paste(unlist(best_hidden_layers), collapse = ", "), "\n")  

nnet_model <- best_model_nnet

#显示模型信息 
summary(nnet_model )


#############4.训练集模型效果评价##############################################

######4.1 模型预测结果####

# Logistic模型
train_prob_lr <- predict(lr_model, newdata = tlog, 
                         type = 'response')     # 指定预测输出为概率,预测Yes的概率
train_prob_lr
train_pred_lr <- factor(ifelse(train_prob_lr > 0.5,'Yes','No'))  # 预测分类
train_pred_lr 

# 决策树模型
train_pred_tree <- predict(tree_model, 
                           newdata = tlog, 
                           type = "class") # 预测分类 
train_pred_tree
train_prob_tree <- predict(tree_model, newdata = tlog, 
                           type = "prob")[, 2]  # 预测Yes的概率
train_prob_tree

# 随机森林
train_pred_rf <- predict(rf_model, newdata = tlog)   # 预测分类
train_pred_rf
train_prob_rf <- predict(rf_model, newdata = tlog, 
                         type = "prob")[, 2]  # 预测Yes的概率  
train_prob_rf

# Xgboost模型
train_prob_xgb <- predict(xgb_model, train_matrix)   # 预测Yes的概率
train_prob_xgb
train_pred_xgb <- factor(ifelse(train_prob_xgb > 0.5,'Yes','No')) # 预测分类
train_pred_xgb

# LightGBM模型
train_prob_lightgbm <- predict(lightgbm_model,
                               newdata = as.matrix(tlog[, selected_vars]),
                               type = 'prob')   # 预测Yes的概率           
train_prob_lightgbm
train_pred_lightgbm <- predict(lightgbm_model,
                               newdata = as.matrix(tlog[, selected_vars]),
                               type = 'class')  # 预测分类
train_pred_lightgbm <- factor(train_pred_lightgbm,levels = c(0,1),labels = c('No','Yes'))
train_pred_lightgbm

# knn 模型
train_pred_knn <- predict(knn_model, newdata = tlog)    # 预测分类
train_pred_knn
train_prob_knn <- predict(knn_model, newdata = tlog, type = "prob")[,"Yes"] # 预测Yes的概率    
train_prob_knn

# 支持向量机
train_pred_svm  <- predict(svm_model, newdata = tlog)   # 预测分类
train_pred_svm
train_prob_svm <- attr(predict(svm_model, newdata = tlog, probability = TRUE), 
                       "probabilities")[, "Yes"]    # 预测Yes的概率
train_prob_svm

# 神经网络 
train_prob_nnet <- predict(nnet_model, tlog)[,2]  # 预测Yes的概率 
train_prob_nnet
train_pred_nnet <- factor(ifelse(train_prob_nnet > 0.5,'Yes','No'))  # 预测分类
train_pred_nnet          


#########4.2 混淆矩阵####

# Logistic模型
confusion_matrix_lr <- caret::confusionMatrix(train_pred_lr, 
                                              tlog$diabetes, 
                                              positive = "Yes")   
print(confusion_matrix_lr) 

# 决策树模型
confusion_matrix_tree <- caret::confusionMatrix(train_pred_tree, 
                                                tlog$diabetes, 
                                                positive = "Yes") # 训练集
print(confusion_matrix_tree) 

# 随机森林
confusion_matrix_rf <- caret::confusionMatrix(train_pred_rf, 
                                              tlog$diabetes,
                                              positive = "Yes")  
print(confusion_matrix_rf) 

# Xgboost模型
confusion_matrix_xgb <- caret::confusionMatrix(train_pred_xgb, 
                                               tlog$diabetes, 
                                               positive = "Yes")  
print(confusion_matrix_xgb) 


# LightGBM模型
confusion_matrix_lightgbm <- caret::confusionMatrix(train_pred_lightgbm, 
                                                    tlog$diabetes, 
                                                    positive = "Yes") 
print(confusion_matrix_lightgbm) 

# knn 模型
confusion_matrix_knn <- caret::confusionMatrix(train_pred_knn, 
                                               tlog$diabetes, 
                                               positive = "Yes") 
print(confusion_matrix_knn) 

# 支持向量机
confusion_matrix_svm <- caret::confusionMatrix(train_pred_svm, 
                                               tlog$diabetes, 
                                               positive = "Yes")
print(confusion_matrix_svm) 

# 神经网络 
confusion_matrix_nnet <- caret::confusionMatrix(train_pred_nnet, 
                                                tlog$diabetes, 
                                                positive = "Yes")  
print(confusion_matrix_nnet)           

#########4.3 ROC曲线####

## 计算ROC的auc值及95%CI
# (1) lr
roc_lr <- roc(tlog$diabetes,    # 目标变量的真实标签
              as.numeric(train_prob_lr)  # 模型预测的概率值
              )
auc_lr <- roc_lr$auc  # AUC值
auc_lr    
ci.auc(roc_lr)        # AUC值的95%CI

# (2) tree
roc_tree <- roc(tlog$diabetes, as.numeric(train_prob_tree))
auc_tree <- roc_tree$auc
auc_tree
ci.auc(auc_tree)

# (3) rf
roc_rf <- roc(tlog$diabetes, as.numeric(train_prob_rf))
auc_rf <- roc_rf$auc
auc_rf
ci.auc(auc_rf)

# (4) xgboost
roc_xgb <- roc(tlog$diabetes, as.numeric(train_prob_xgb))
auc_xgb <- roc_xgb$auc
auc_xgb
ci.auc(auc_xgb)

# (5) lightgbm
roc_lightgbm <- roc(tlog$diabetes, as.numeric(train_prob_lightgbm))
auc_lightgbm <- roc_lightgbm$auc
auc_lightgbm
ci.auc(auc_lightgbm)

# (6) knn
roc_knn <- roc(tlog$diabetes, as.numeric(train_prob_knn))
auc_knn <- roc_knn$auc
auc_knn
ci.auc(auc_knn)

# (7) svm
roc_svm <- roc(tlog$diabetes, as.numeric(train_prob_svm))
auc_svm <- roc_svm$auc
auc_svm
ci.auc(auc_svm)

# (8) nnet
roc_nnet <- roc(tlog$diabetes, as.numeric(train_prob_nnet))
auc_nnet <- roc_nnet$auc
auc_nnet
ci.auc(auc_nnet)

# 绘制ROC曲线
plot(roc_lr, 
     col = "red", # 曲线颜色为红色
     lwd = 2,     # 曲线的线宽为 2
     main = "ROC Curves for Training dataset",  # 设置图的标题
     xlab = "1 - Specificity", ylab = "Sensitivity", # 设置X轴和Y轴标签
     legacy.axes = TRUE,  # 使 X 轴范围从 0 到 1
     cex.main = 1.5,      # 设置标题字体大小
     cex.lab = 1.2, cex.axis = 1.2  # 设置坐标轴标签和刻度字体大小
     )
lines(roc_tree, col = "blue", lwd = 2)
lines(roc_rf, col = "green", lwd = 2)
lines(roc_knn, col = "purple", lwd = 2)
lines(roc_svm, col = "orange", lwd = 2)
lines(roc_nnet, col = "brown", lwd = 2)
lines(roc_xgb, col = "pink", lwd = 2)
lines(roc_lightgbm, col = "cyan", lwd = 2)

# 添加图例
legend("bottomright", 
       legend = c(    # 创建每个模型的AUC文本
         paste("Logistic Regression (AUC = ", round(auc_lr, 3), ")", sep = ""),
         paste("Decision Tree (AUC = ", round(auc_tree, 3), ")", sep = ""),
         paste("Random Forest (AUC = ", round(auc_rf, 3), ")", sep = ""),
         paste("KNN (AUC = ", round(auc_knn, 3), ")", sep = ""),
         paste("SVM (AUC = ", round(auc_svm, 3), ")", sep = ""),
         paste("Neural Network (AUC = ", round(auc_nnet, 3), ")", sep = ""),
         paste("XGBoost (AUC = ", round(auc_xgb, 3), ")", sep = ""),
         paste("LightGBM (AUC = ", round(auc_lightgbm, 3), ")", sep = "") ), 
       col = c("red", "blue", "green", "purple", "orange", 
               "brown", "pink", "cyan"), 
       lty = 1,   # 设置图例中线的样式为实线
       lwd = 2,   # 设置图例中线的宽度为 2
       cex = 0.6) # 设置图例字体的缩放比例


#############4.4 校准曲线####

# 预测结果及真实标签汇总为一个数据框
calibration_data <- data.frame(
  Model = c(rep("Logistic Regression", length(train_prob_lr)),# rep():重复指定字符串多次
            rep("Decision Tree", length(train_prob_tree)),
            rep("Random Forest", length(train_prob_rf)),
            rep("KNN", length(train_prob_knn)),
            rep("SVM", length(train_prob_svm)),
            rep("Neural Network", length(train_prob_nnet)),
            rep("XGBoost", length(train_prob_xgb)),
            rep("LightGBM", length(train_prob_lightgbm))),
  Probability = c(train_prob_lr, 
                  train_prob_tree,
                  train_prob_rf,
                  train_prob_knn,
                  train_prob_svm,
                  train_prob_nnet,
                  train_prob_xgb,
                  train_prob_lightgbm),#将所有模型的预测概率按顺序拼接成一个向量
  Actual = as.numeric(c(tlog$diabetes)) - 1  # 将因子转为数值
  )

# 绘制光滑校准曲线
ggplot(calibration_data, aes(x = Probability, y = Actual, color = Model)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.5) +  # 使用 LOESS 光滑曲线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +  # 理想参考线
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "Calibration Curves",
    x = "Actual Probability",
    y = "Observed Proportion" ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(colour = "black")) +
  scale_color_brewer(palette = "Set1")  # 使用醒目的配色方案

############4.5 决策曲线分析(DCA)曲线######

# 将模型的预测结果和真实标签组合成一个数据框
dca_tlog <- data.frame(diabetes = as.numeric(tlog$diabetes)-1, 
                       train_prob_lr,
                       train_prob_tree, 
                       train_prob_rf,
                       train_prob_knn,
                       train_prob_svm,
                       train_prob_nnet,
                       train_prob_xgb,
                       train_prob_lightgbm)

# Logistic 模型
dca.result_lr <- decision_curve(diabetes ~ train_prob_lr, 
                                data = dca_tlog, 
                                bootstraps = 10)
# 决策树
dca.result_tree <- decision_curve(diabetes ~ train_prob_tree, 
                                  data = dca_tlog,
                                  bootstraps = 10)
# 随机森林
dca.result_rf <- decision_curve(diabetes ~ train_prob_rf, 
                                data = dca_tlog, 
                                bootstraps = 10)
# XGBoost 
dca.result_xgb <- decision_curve(diabetes ~ train_prob_xgb, 
                                 data = dca_tlog, 
                                 bootstraps = 10)
# LightGBM
dca.result_lightgbm <- decision_curve(diabetes ~ train_prob_lightgbm, 
                                      data = dca_tlog, 
                                      bootstraps = 10)
# knn
dca.result_knn <- decision_curve(diabetes ~ train_prob_knn, 
                                 data = dca_tlog,
                                 bootstraps = 10)
# 支持向量机(SVM)
dca.result_svm <- decision_curve(diabetes ~ train_prob_svm, 
                                 data = dca_tlog)
# 神经网络
dca.result_nnet <- decision_curve(diabetes ~ train_prob_nnet, 
                                  data = dca_tlog, 
                                  bootstraps = 10)

plot_decision_curve(
  list(dca.result_lr, dca.result_tree, dca.result_rf, 
       dca.result_knn, dca.result_svm, dca.result_nnet,
       dca.result_xgb, dca.result_lightgbm),  # 传入包含决策曲线分析结果的列表
  curve.names = c("Logistic Regression", "Decision Tree", "Random Forest", 
                  "KNN", "SVM", "Neural Network", "XGBoost", "LightGBM"),
  col = c("red", "green", "blue", "purple", "orange", "brown", "cyan", "magenta"),
  lwd = 2,  # 设置线宽
  confidence.intervals = FALSE, # 禁用置信区间
  legend.position = ("none")  
)

# 手动添加图例
legend("bottomleft", 
       legend = c("Logistic Regression", "Decision Tree", "Random Forest", 
                  "KNN", "SVM", "Neural Network", "XGBoost", "LightGBM"),
       col = c("red", "green", "blue", "purple", "orange", "brown", "cyan", "magenta"),
       lwd = 2,
       cex = 0.6,  # 调整cex改变字体大小
       bty = "y",   # 无边框
       y.intersp = 0.8,    # 调整条目之间的间距
       x.intersp = 0.5,    # 调整颜色线条与标签之间的水平间距
       text.width = 0.1)  # 自动调整宽度
       

##############################5.验证集模型效果评价##########################

######5.1 模型预测结果####

# 结局变量因子化
valdata$diabetes <- factor(valdata$diabetes,levels = c(0,1),labels = c('No','Yes'))

# Logistic模型
val_prob_lr <- predict(lr_model, newdata = valdata, 
                        type = 'response')    # 预测Yes的概率
val_prob_lr
val_pred_lr <-  factor(ifelse(val_prob_lr > 0.5,'Yes','No'))  # 预测分类
val_pred_lr 

# 决策树模型
val_pred_tree <- predict(tree_model, newdata = valdata, type = "class") # 预测分类 
val_pred_tree
val_prob_tree <- predict(tree_model, newdata = valdata, 
                           type = "prob")[, 2]  # 预测Yes的概率
val_prob_tree


# 随机森林
val_pred_rf <- predict(rf_model, newdata = valdata)   # 预测分类
val_pred_rf
val_prob_rf <- predict(rf_model, newdata = valdata, 
                         type = "prob")[, 2]  # 预测Yes的概率  
val_prob_rf

# Xgboost模型
val_prob_xgb <- predict(xgb_model, val_matrix)   # 预测Yes的概率
val_prob_xgb
val_pred_xgb <- factor(ifelse(val_prob_xgb > 0.5,'Yes','No')) # 预测分类
val_pred_xgb

# LightGBM模型
val_prob_lightgbm <- predict(lightgbm_model,
                               newdata = as.matrix(valdata[,selected_vars]),
                               type = 'prob')   # 预测Yes的概率           
val_prob_lightgbm
val_pred_lightgbm <- predict(lightgbm_model,
                               newdata = as.matrix(valdata[,selected_vars]),
                               type = 'class')  # 预测分类
val_pred_lightgbm <- factor(val_pred_lightgbm,
                             levels = c(0,1),labels = c('No','Yes'))
val_pred_lightgbm

# knn 模型
val_pred_knn <- predict(knn_model, newdata = valdata) # 预测分类
val_pred_knn
val_prob_knn <- predict(knn_model, newdata = valdata, 
                         type = "prob")[,"Yes"]  # 预测Yes的概率    
val_prob_knn

# 支持向量机
val_pred_svm  <- predict(svm_model, newdata = valdata)   # 预测分类
val_pred_svm
val_prob_svm <- attr(predict(svm_model, newdata = valdata, probability = TRUE), 
                       "probabilities")[, "Yes"]    # 预测Yes的概率
val_prob_svm

# 神经网络 
val_prob_nnet <- predict(nnet_model, valdata)[,2]  # 预测Yes的概率 
val_prob_nnet
val_pred_nnet <- factor(ifelse(val_prob_nnet > 0.5,'Yes','No'))  # 预测分类
val_pred_nnet          

#########5.2 混淆矩阵####

# Logistic模型
confusion_matrix_lr1 <- caret::confusionMatrix(
  val_pred_lr, valdata$diabetes, positive = "Yes")   
print(confusion_matrix_lr1) 

# 决策树模型
confusion_matrix_tree1 <- caret::confusionMatrix(val_pred_tree, 
                                                valdata$diabetes, 
                                                positive = "Yes") # 训练集
print(confusion_matrix_tree1) 

# 随机森林
confusion_matrix_rf1 <- caret::confusionMatrix(val_pred_rf, 
                                              valdata$diabetes,
                                              positive = "Yes")  
print(confusion_matrix_rf1) 

# Xgboost模型
confusion_matrix_xgb1 <- caret::confusionMatrix(val_pred_xgb, 
                                               valdata$diabetes, 
                                               positive = "Yes")  
print(confusion_matrix_xgb1) 


# LightGBM模型
confusion_matrix_lightgbm1 <- caret::confusionMatrix(val_pred_lightgbm, 
                                                    valdata$diabetes, 
                                                    positive = "Yes") 
print(confusion_matrix_lightgbm1) 

# knn 模型
confusion_matrix_knn1 <- caret::confusionMatrix(val_pred_knn, 
                                               valdata$diabetes, 
                                               positive = "Yes") 
print(confusion_matrix_knn1) 

# 支持向量机
confusion_matrix_svm1 <- caret::confusionMatrix(val_pred_svm, 
                                               valdata$diabetes, 
                                               positive = "Yes")
print(confusion_matrix_svm1) 

# 神经网络 
confusion_matrix_nnet1 <- caret::confusionMatrix(val_pred_nnet, 
                                                valdata$diabetes, 
                                                positive = "Yes")  
print(confusion_matrix_nnet1)

#########5.3 ROC曲线####

# 计算ROC的auc值及95%CI

roc_lr_val <- roc(valdata$diabetes, as.numeric(val_prob_lr))
auc_lr_val <- roc_lr_val$auc  # AUC值
auc_lr_val    
ci.auc(roc_lr_val)        # AUC值的95%CI

roc_tree_val <- roc(valdata$diabetes, as.numeric(val_prob_tree))
auc_tree_val <- roc_tree_val $auc
auc_tree_val
ci.auc(auc_tree_val)

roc_rf_val <- roc(valdata$diabetes, as.numeric(val_prob_rf))
auc_rf_val <- roc_rf_val $auc
auc_rf_val
ci.auc(auc_rf_val)

roc_xgb_val <- roc(valdata$diabetes, as.numeric(val_prob_xgb))
auc_xgb_val <- roc_xgb_val $auc
auc_xgb_val
ci.auc(auc_xgb_val)

roc_lightgbm_val <- roc(valdata$diabetes, as.numeric(val_prob_lightgbm))
auc_lightgbm_val <- roc_lightgbm_val $auc
auc_lightgbm_val
ci.auc(auc_lightgbm_val)

roc_knn_val <- roc(valdata$diabetes, as.numeric(val_prob_knn))
auc_knn_val <- roc_knn_val $auc
auc_knn_val
ci.auc(auc_knn_val)

roc_svm_val <- roc(valdata$diabetes, as.numeric(val_prob_svm))
auc_svm_val <- roc_svm_val$auc
auc_svm_val
ci.auc(auc_svm_val)

roc_nnet_val <- roc(valdata$diabetes, as.numeric(val_prob_nnet))
auc_nnet_val <- roc_nnet_val $auc
auc_nnet_val
ci.auc(auc_nnet_val)

# 绘制ROC曲线
plot(roc_lr_val, col = "red", lwd = 2, main = "ROC Curves for val dataset", 
     xlab = "1 - Specificity", ylab = "Sensitivity", legacy.axes = TRUE, 
     cex.main = 1.6, cex.lab = 1.3, cex.axis = 1.2)
lines(roc_tree_val, col = "blue", lwd = 2)
lines(roc_rf_val, col = "green", lwd = 2)
lines(roc_knn_val, col = "purple", lwd = 2)
lines(roc_svm_val, col = "orange", lwd = 2)
lines(roc_nnet_val, col = "brown", lwd = 2)
lines(roc_xgb_val, col = "pink", lwd = 2)
lines(roc_lightgbm_val, col = "cyan", lwd = 2)

# 添加图例
legend("bottomright", 
       legend = c(
         paste("Logistic Regression (AUC = ", round(auc_lr_val, 3), ")", sep = ""),
         paste("Decision Tree (AUC = ", round(auc_tree_val, 3), ")", sep = ""),
         paste("Random Forest (AUC = ", round(auc_rf_val, 3), ")", sep = ""),
         paste("KNN (AUC = ", round(auc_knn_val, 3), ")", sep = ""),
         paste("SVM (AUC = ", round(auc_svm_val, 3), ")", sep = ""),
         paste("Neural Network (AUC = ", round(auc_nnet_val, 3), ")", sep = ""),
         paste("XGBoost (AUC = ", round(auc_xgb_val, 3), ")", sep = ""),
         paste("LightGBM (AUC = ", round(auc_lightgbm_val, 3), ")", sep = "") ), 
       col = c("red", "blue", "green", "purple", "orange", 
               "brown", "pink", "cyan"), 
       lty = 1, lwd = 2, cex = 0.6)


#############5.4 校准曲线####

# 预测结果及真实标签汇总为一个数据框
calibration_data1 <- data.frame(
  Model = c(rep("Logistic Regression", length(val_prob_lr)),
            rep("Decision Tree", length(val_prob_tree)),
            rep("Random Forest", length(val_prob_rf)),
            rep("KNN", length(val_prob_knn)),
            rep("SVM", length(val_prob_svm)),
            rep("Neural Network", length(val_prob_nnet)),
            rep("XGBoost", length(val_prob_xgb)),
            rep("LightGBM", length(val_prob_lightgbm))),
  Probability = c(val_prob_lr, 
                  val_prob_tree,
                  val_prob_rf,
                  val_prob_knn,
                  val_prob_svm,
                  val_prob_nnet,
                  val_prob_xgb,
                  val_prob_lightgbm),
  Actual = as.numeric(c(valdata$diabetes)) - 1 )  # 将因子转为数值
 


# 绘制光滑校准曲线
ggplot(calibration_data1, aes(x = Probability, y = Actual, color = Model)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.5) +  # 使用 LOESS 光滑曲线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +  # 理想参考线
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = "Calibration Curves",
    x = "Actual Probability",
    y = "Observed Proportion" ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(colour = "black")) +
  scale_color_brewer(palette = "Set1")  # 使用醒目的配色方案

############5.5 DCA曲线######

dca_valdata <- data.frame(diabetes = as.numeric(valdata$diabetes)-1, 
                       val_prob_lr,
                       val_prob_tree, 
                       val_prob_rf,
                       val_prob_knn,
                       val_prob_svm,
                       val_prob_nnet,
                       val_prob_xgb,
                       val_prob_lightgbm)
# Logistic Regression 
dca.result_lr1 <- decision_curve(diabetes ~ val_prob_lr, 
                                data = dca_valdata, 
                                bootstraps = 10)
# Decision Tree
dca.result_tree1 <- decision_curve(diabetes ~ val_prob_tree, 
                                  data = dca_valdata, 
                                  bootstraps = 10)
# Random Forest
dca.result_rf1 <- decision_curve(diabetes ~ val_prob_rf, 
                                data = dca_valdata, 
                                bootstraps = 10)
# XGBoost
dca.result_xgb1 <- decision_curve(diabetes ~ val_prob_xgb, 
                                 data = dca_valdata, 
                                 bootstraps = 10)
# LightGBM
dca.result_lightgbm1 <- decision_curve(diabetes ~ val_prob_lightgbm, 
                                      data = dca_valdata, 
                                      bootstraps = 10)
# KNN
dca.result_knn1 <- decision_curve(diabetes ~ val_prob_knn, 
                                 data = dca_valdata, 
                                 bootstraps = 10)
# SVM
dca.result_svm1 <- decision_curve(diabetes ~ val_prob_svm, 
                                 data = dca_valdata, 
                                 bootstraps = 10)
# Neural Network
dca.result_nnet1 <- decision_curve(diabetes ~ val_prob_nnet, 
                                  data = dca_valdata, 
                                  bootstraps = 10)

plot_decision_curve(
  list(dca.result_lr1, dca.result_tree1, dca.result_rf1, 
       dca.result_knn1, dca.result_svm1, dca.result_nnet1,
       dca.result_xgb1, dca.result_lightgbm1),  # 传入包含决策曲线分析结果的列表
  curve.names = c("Logistic Regression", "Decision Tree", "Random Forest", 
                  "KNN", "SVM", "Neural Network", "XGBoost", "LightGBM"),
  col = c("red", "green", "blue", "purple", "orange", "brown", "cyan", "magenta"),
  lwd = 2,  # 设置线宽
  confidence.intervals = FALSE, # 禁用置信区间
  legend.position = ("none")  
)

# 手动添加图例
legend("bottomright", 
       legend = c("Logistic Regression", "Decision Tree", "Random Forest", 
                  "KNN", "SVM", "Neural Network", "XGBoost", "LightGBM"),
       col = c("red", "green", "blue", "purple", "orange", "brown", "cyan", "magenta"),
       lwd = 2,
       cex = 0.6,  # 通过调整cex改变字体大小
       bty = "y",   # 无边框
       y.intersp = 0.8,    # 调整条目之间的间距
       x.intersp = 0.5,    # 调整颜色线条与标签之间的水平间距
       text.width = 0.15)  # 自动调整宽度


#####################6.模型的SHAP解释#########################################

######6.1 xgboost#### 

# shap解释使用数据量(按照自己需要求更改,不能超过验证集的样本量)
n_tlog = 150    
n_valdata = 150  

# 数据转换为矩阵
tlog_matrix <- as.matrix(tlog[1:n_tlog, selected_vars])  # 训练集中提取自变量  
valdata_matrix <- as.matrix(valdata[1:n_valdata, selected_vars])  # 验证集中提取自变量  

# 使用 kernelshap 函数来计算 SHAP 值  
explain_kernel_xgb <- kernelshap(xgb_model,     # 训练好的xgboost模型
                                 tlog_matrix,   # 训练集数据
                                 bg_X = valdata_matrix)  # 使用验证集的自变量
# 通过 shapviz 生成 SHAP 值的可视化对象
shap_value_xgb <- shapviz(explain_kernel_xgb,
                      X_pred = tlog[1:n_tlog,selected_vars],   # 用训练集自变量生成 SHAP 可视化
                      interactions = TRUE)     # 启用特征间的交互效应

# (1) 特征重要性条形图
sv_importance(shap_value_xgb,     # SHAP 值对象
              kind = "bar",       # 指定图的类型为条形图
              show_numbers = T,   # 不显示特征重要性值
              fill = "#1f77b4"     # 设置条形图的颜色为 #1f77b4（蓝色）
) + 
  theme_bw() +  # 设置黑白主题
  ggtitle("xgboost") +  # 设置标题为 xgboost
  theme(plot.title = element_text(hjust = 0.5,  # 使标题居中
                                  face = "bold",  # 标题字体加粗
                                  color = "black"))  # 标题字体颜色为黑色

# (2) 特征重要性蜂群图
sv_importance(shap_value_xgb, 
              kind = "beeswarm")

sv_importance(shap_value_xgb, 
              kind = "beeswarm",   # 指定图的类型为蜂群图
              viridis_args = list(begin =0.2, end =0.9, option ="B"), # 设置颜色方案
              show_numbers = F)+  # 不显示数字
  ggtitle("xgboost")+   # 标题为 xgboost
  theme_bw()+           # 黑白背景主题
  theme(plot.title = element_text(hjust = 0.5, face ="bold", color ="black")) 

# 两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_xgb, kind = "both") + theme_bw()


# (3) 依赖图

# 分析glucose变量与SHAP值的关系,并考虑 pregnant 变量
sv_dependence(shap_value_xgb,
              v = "glucose",  #指定第一个变量,横坐标
              color_var = "pregnant",    # 第二个变量,区分点的颜色,展示交互信息
              size = 4)+                 # 点的大小设置为 4
  theme_bw()+
  ggtitle("xgboost")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# 分析age变量与SHAP值的关系,并考虑 exercise 变量

# exercise是数值型的,先将它转换为因子类型
shap_value_xgb$X$exercise <- factor(shap_value_xgb$X$exercise, levels = c(0, 1))

sv_dependence(shap_value_xgb,
              v = "pregnant",  #指定第一个变量
              color_var = "exercise",    #指定第二个变量
              size = 4) + 
  theme_bw()+
  ggtitle("xgboost")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black")) +
  scale_color_manual(values = c("0" = "#1f77b4", "1" = "#ff6347"),  # 设置颜色
                     labels = c("No", "Yes"))  # 设置图例标签

# (4) 单样本特征：瀑布图
sv_waterfall(shap_value_xgb, 
             row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("xgboost")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征：力图
sv_force(shap_value_xgb, row_id = 12,size = 9)+
  ggtitle("xgboost")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#############6.2 RF模型#######

explain_kernel_rf <- kernelshap(rf_model, 
                                tlog[1:n_tlog,selected_vars], 
                                bg_X = valdata[1:n_valdata,selected_vars],
                                pred_fun = function(model, X)
                                  predict(model, X, type = "prob")[, 2])  
shap_value_rf <- shapviz(explain_kernel_rf, 
                         X_pred = tlog[1:n_tlog,selected_vars], 
                         interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_rf, 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4")+
  theme_bw()+
  ggtitle("randomForest") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (2) 特征重要性蜂群图
sv_importance(shap_value_rf, 
              kind = "beeswarm", 
              viridis_args = list(begin=0.25, end=0.85, option="B"),
              show_numbers = F)+
  ggtitle("randomForest")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_rf, kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_rf,
              v = "glucose", # 指定第一个变量
              color_var = "pregnant",  # 指定第二个变量
              size = 4)+
  theme_bw()+
  ggtitle("randomForest")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_rf, row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("randomForest")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_rf, row_id = 12,size = 9) + 
  ggtitle("randomForest")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))


######6.3 logistic模型####

explain_kernel_lr <- kernelshap(lr_model, 
                                tlog[1:n_tlog,selected_vars], 
                                bg_X = valdata[1:n_valdata,selected_vars])  
shap_value_lr <- shapviz(explain_kernel_lr,
                         X_pred = tlog[1:n_tlog,selected_vars], 
                         interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_lr, 
              kind = "bar", 
              show_numbers = T,  # T显示特征重要性值,F 不显示
              fill = "#1f77b4" )+
  theme_bw()+
  ggtitle("Logistic")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))


# (2) 特征重要性蜂群图
sv_importance(shap_value_lr, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),
              show_numbers = FALSE)+
  ggtitle("Logistic")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_lr, kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_lr,
              v = "glucose",        #指定第一个变量
              color_var = "pregnant", #指定第二个变量
              size = 4)+
  theme_bw()+
  ggtitle("Logistic")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_lr,
             row_id = 2,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("Logistic")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_lr, row_id = 2,size = 10)+
  ggtitle(label = "Logistic")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))


######6.4 决策树模型####

explain_kernel_tree <- kernelshap(tree_model, 
                                  tlog[1:n_tlog,selected_vars], 
                                  bg_X = valdata[1:n_valdata,selected_vars])  
shap_value_tree <- shapviz(explain_kernel_tree,
                           X_pred = tlog[1:n_tlog,selected_vars], 
                           interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_tree$Yes, 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4",
              class = "Yes")+
  theme_bw()+
  ggtitle("decision tree")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (2) 特征重要性蜂群图
sv_importance(shap_value_tree$Yes, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),#A-H
              show_numbers = F)+
  ggtitle("decision tree")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_tree$Yes, kind = "both") + theme_bw()


# (3) 依赖图
sv_dependence(shap_value_tree$Yes,
              v = "glucose",       #指定第1个变量
              color_var = "pregnant",  #指定第2个变量
              size = 4)+
  theme_bw()+
  ggtitle("decision tree")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))


# (4) 瀑布图
sv_waterfall(shap_value_tree$Yes, 
             row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("decision tree")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_tree$Yes, row_id = 12,size = 9)+
  ggtitle("decision tree")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#####6.5 knn模型####

explain_kernel_knn <- kernelshap(knn_model, 
                                 tlog[1:n_valdata,selected_vars_scaled], 
                                 bg_X = valdata[1:n_valdata,selected_vars_scaled],
                                 pred_fun = function(model, X)
                                   predict(model, X, type = "prob")[, "Yes"])  
shap_value_knn <- shapviz(explain_kernel_knn,
                          X_pred = tlog_scaled[1:n_tlog,selected_vars_scaled], 
                          interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_knn, 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4")+
  theme_bw()+
  ggtitle("knn")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (2) 特征重要性蜂群图
sv_importance(shap_value_knn, kind = "beeswarm")

sv_importance(shap_value_knn, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),#A-H
              show_numbers = F)+
  ggtitle("knn")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_knn, kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_knn,v = "glucose_scaled",  #指定第一个变量
              color_var = "pregnant",   #指定第二个变量
              size = 4)+
  theme_bw()+
  ggtitle("knn")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_knn, 
             row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("knn")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_knn, row_id = 12,size = 9)+
  ggtitle("knn")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#####6.6 LightGBM模型####

tlog_matrix <- as.matrix(tlog[1:n_tlog, selected_vars])  # 去掉结果列  
bg_X_matrix <- as.matrix(valdata[1:n_valdata, selected_vars])  # 同样处理验证集  

# 使用 kernelshap 进行解释  
explain_kernel_lightgbm <- kernelshap(lightgbm_model, 
                                      tlog_matrix, 
                                      bg_X = bg_X_matrix)  
shap_value_lightgbm <- shapviz(explain_kernel_lightgbm,
                               X_pred = tlog[1:n_tlog,selected_vars], 
                               interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_lightgbm, 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4" )+
  theme_bw()+
  ggtitle("LightGBM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (2) 特征重要性蜂群图
sv_importance(shap_value_lightgbm, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),
              show_numbers = F)+
  ggtitle("LightGBM")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# 两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_lightgbm, kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_lightgbm, 
              v = "glucose",  #指定第一个变量
              color_var = "pregnant",  #指定第二个变量
              size = 4 )+
  theme_bw()+
  ggtitle("LightGBM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_lightgbm, row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("LightGBM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_lightgbm, row_id = 12,size = 9)+
  ggtitle("LightGBM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

######6.7 SVM模型####

explain_kernel_svm <- kernelshap(svm_model, 
                                 tlog[1:n_tlog,selected_vars_scaled], 
                                 bg_X = valdata[1:n_valdata,selected_vars_scaled],
                                 pred_fun = function(model, X)
              attr(predict(model, X, probability=TRUE),"probabilities")[,"Yes"]) 

shap_value_svm  <- shapviz(explain_kernel_svm ,
                           X_pred = tlog[1:n_tlog,selected_vars_scaled], 
                           interactions = TRUE) 


# (1) 特征重要性条形图
sv_importance(shap_value_svm , 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4")+
  theme_bw()+
  ggtitle("SVM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (2) 特征重要性蜂群图
sv_importance(shap_value_svm, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "C"),
              show_numbers = F)+
  ggtitle("SVM")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_svm , kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_svm ,
              v = "glucose_scaled",   #指定第一个变量
              color_var = "pregnant",    #指定第二个变量
              size = 4)+
  theme_bw()+
  ggtitle("SVM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_svm , row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("SVM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_svm , row_id = 12,size = 9)+
  ggtitle("SVM")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))


#####6.8 nnet模型####

explain_kernel_nnet <- kernelshap(nnet_model, 
                                  tlog[1:n_tlog,selected_vars_scaled], 
                                  bg_X = valdata[1:n_valdata,selected_vars_scaled])  
shap_value_nnet <- shapviz(explain_kernel_nnet ,
                           X_pred = tlog[1:n_tlog, selected_vars_scaled], 
                           interactions = TRUE) 

# (1) 特征重要性条形图
sv_importance(shap_value_nnet$Class_2, 
              kind = "bar", 
              show_numbers = F,
              fill = "#1f77b4" )+
  theme_bw()+
  ggtitle("nnet")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))             

# (2) 特征重要性蜂群图
sv_importance(shap_value_nnet$Class_2, 
              kind = "beeswarm", 
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),
              show_numbers = F)+
  ggtitle("nnet")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

#两种图形(条形图与蜂群图)叠加
sv_importance(shap_value_nnet$Class_2, kind = "both") + theme_bw()

# (3) 依赖图
sv_dependence(shap_value_nnet$Class_2,v = "glucose_scaled", #指定第一个变量
              color_var = "pregnant",  #指定第二个变量
              size = 4)+
  theme_bw()+
  ggtitle("nnet")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (4) 瀑布图
sv_waterfall(shap_value_nnet$Class_2, row_id = 12,
             fill_colors = c("#f7d13d", "#a52c60"))+
  theme_bw()+
  ggtitle("nnet")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))

# (5) 单样本特征
sv_force(shap_value_nnet$Class_2, row_id = 12,size = 9)+
  ggtitle("nnet")+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",color = "black"))



