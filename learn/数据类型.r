# 向量
x <- c(1,2,3,4,5,6)
x <- c("one", "two", "three", "four", "five", "six")

# 矩阵
# 生成5x4的矩阵，元素必须是同种类型，nrow指定行数，ncol指定列数
x <- matrix(1:20, nrow=5, ncol=4)

x <- c("R1", "R2", "R3", "R4", "R5")
y <- c("C1", "C2", "C3", "C4")
z <- matrix(1:20, nrow=5, ncol=4, byrow=TRUE, dimnames=list(x, y))
# byrow=TRUE表示按行填充，FALSE按列填充
# dimnames=list(x, y)，行和列分别按照x和y进行命名
z[1, c(4,5)] # 输出第一行4，5位置的元素

# 数组
# 与矩阵类似，不过维度可以大于2
dim1 <- c("A1", "A2")
dim2 <- c("B1", "B2", "B3")
dim3 <- c("C1", "C2", "C3", "C4")
x <- array(1:24, c(2, 3, 4), dimnames=list(dim1, dim2, dim3))
# 以向量dim1, dim2, dim3作为三个坐标轴，数组的大小为2x3x4,即有2行、3列和4深度，数组的内容是1到24的整数

# 数据框
patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type3", "Type4")
status <- c("Poor", "Improved", "Excellent", "Poor")
patientdata <- data.frame(patientID, age, diabetes, status)
# 每一列数据的类型是唯一的

## 数据框选取数据
patientdata$age
# 使用with()简化提取数据，括号内的赋值仅在括号内生效，使用 <<- 可以使括号内赋值在括号外生效
with(patientdata, {
    summary(age) # 计算变量的汇总统计信息
    plot(patientID, age)
})

# 因子（！！！！未掌握）
x <- c("type1", "type2", "type1", "type1")
x <- factor(x)
# 将以上向量存储为(1,2,1,1)
# factor函数以整数向量的形式存储类别值，整数的取值范围是1-k，k是名义变量中唯一值的个数
# 字符型变量，因子的水平会依照字母顺序创建

# 列表
g <- "My list"
h <- c(23,32,45,43)
j <- matrix(1:10,nrow=5)
k <- c("one", "two", "three")
mylist <- list(title=g,ages=h,j,k)
mylist[[2]] # 访问指定列表

# tibble数据框
mtcars <- as_tibble(mtcars)
# tibble数据框不会将字符变量转换为因子，不会改变变量的名称（R中的命名不能有空格）
# tibble数据框不支持行名