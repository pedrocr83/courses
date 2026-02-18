# Programming Primer for Bioinformatics

**Estimated Time:** 10-20 hours (self-paced)  
**Goal:** Build programming skills to pass the Programming Assessment (7/10)

---

## How to Use This Primer

- Work through sections relevant to your weak areas
- Run all code examples yourself
- Complete practice exercises at the end of each section

---

## 1. Basic Syntax & Operator Precedence

### Python

Variables store values. Assignment uses `=`. Arithmetic follows **PEMDAS** (Parentheses, Exponents, Multiplication/Division, Addition/Subtraction).

```python
# Variables and assignment
gene_count = 25000
expression_value = 5.3

# PEMDAS: multiplication happens before addition
# 10 * 2 = 20, then 5 + 20 = 25
x = 5
y = 10
print(x + y * 2)  # Prints 25, not 30

# Use parentheses to change order
print((x + y) * 2)  # Prints 30
```

```python
# Bioinformatics example: calculate fold-change
baseline = 10
treated = 50
fold_change = treated / baseline  # 5.0
print(f"Fold change: {fold_change}x")
```

### R

R uses `<-` for assignment (though `=` also works). Same PEMDAS rules apply.

```r
# Variables and assignment
gene_count <- 25000
expression_value <- 5.3

# PEMDAS: multiplication before addition
x <- 5
y <- 10
print(x + y * 2)  # Prints 25, not 30

# Parentheses change order
print((x + y) * 2)  # Prints 30
```

```r
# Bioinformatics example: fold-change
baseline <- 10
treated <- 50
fold_change <- treated / baseline
print(paste("Fold change:", fold_change, "x"))
```

**Practice:** Predict the output of `print(2 + 3 * 4)` and `print((2 + 3) * 4)` in both languages.

---

## 2. Data Structures: Lists (Python) / Vectors (R)

### Python: Lists

Lists hold ordered collections. Indexing starts at 0. Use slicing for subsequences.

```python
# Create list of expression values (TPM for 5 genes)
expression = [10.2, 5.1, 0.3, 22.0, 8.7]
print(expression[0])      # First: 10.2
print(expression[-1])     # Last: 8.7
print(expression[1:4])    # Slice: [5.1, 0.3, 22.0]

# Iterate and filter even numbers
numbers = list(range(1, 11))
for num in numbers:
    if num % 2 == 0:
        print(num)
```

```python
# List comprehension for even numbers (compact)
evens = [n for n in range(1, 11) if n % 2 == 0]
print(evens)  # [2, 4, 6, 8, 10]
```

### R: Vectors

Vectors are the fundamental R structure. Indexing starts at 1.

```r
# Create vector of expression values
expression <- c(10.2, 5.1, 0.3, 22.0, 8.7)
print(expression[1])      # First: 10.2
print(expression[length(expression)])  # Last: 8.7
print(expression[2:4])    # Elements 2-4

# Filter even numbers from 1 to 10
numbers <- 1:10
evens <- numbers[numbers %% 2 == 0]
print(evens)  # 2 4 6 8 10
```

```r
# Using a for loop
for (num in 1:10) {
  if (num %% 2 == 0) {
    print(num)
  }
}
```

**Practice:** Create a list/vector of numbers 1–10 and print only odd numbers.

---

## 3. Reading Data Files

### Python: pandas

```python
import pandas as pd

# Read CSV (e.g., gene expression matrix)
df = pd.read_csv("expression_data.csv")

# Inspect the data
print(df.head())           # First 5 rows
print(df.shape)            # (rows, columns)
print(df.describe())       # Summary statistics
```

```python
# If file has different delimiter or no header
# df = pd.read_csv("data.txt", sep="\t", header=0)
```

### R: read.csv

```r
# Read CSV
df <- read.csv("expression_data.csv")

# Inspect the data
head(df)                   # First 6 rows
dim(df)                    # c(rows, columns)
summary(df)                # Summary statistics
```

```r
# Alternative: read.table for tab-separated
# df <- read.table("data.txt", sep="\t", header=TRUE)
```

**Practice:** Create a small CSV with columns `gene`, `expression`, `pvalue`. Load it and print its dimensions.

---

## 4. Data Manipulation

### Python: Filtering and Selecting

```python
import pandas as pd

# Filter rows where score > 50
filtered = df[df['score'] > 50]

# Select specific columns
genes_only = df[['gene', 'expression']]

# Sort by expression (descending)
sorted_df = df.sort_values('expression', ascending=False)

# Group by condition and get mean
grouped = df.groupby('condition')['expression'].mean()
```

### R: Filtering and Selecting

```r
# Filter rows where score > 50
filtered <- df[df$score > 50, ]

# Select columns
genes_only <- df[, c("gene", "expression")]

# Sort by expression
sorted_df <- df[order(df$expression, decreasing = TRUE), ]

# With dplyr (if installed):
# library(dplyr)
# filtered <- filter(df, score > 50)
# grouped <- df %>% group_by(condition) %>% summarise(mean_exp = mean(expression))
```

**Practice:** Given a DataFrame with `pvalue`, filter rows where `pvalue < 0.05`.

---

## 5. Plotting Basics

### Python: matplotlib and seaborn

```python
import matplotlib.pyplot as plt

# Scatter plot: expression vs p-value
plt.scatter(df['expression'], df['pvalue'])
plt.xlabel('Expression (TPM)')
plt.ylabel('P-value')
plt.title('Gene Expression vs Significance')
plt.show()
```

```python
# Histogram of expression values
plt.hist(df['expression'], bins=20, edgecolor='black')
plt.xlabel('Expression')
plt.ylabel('Count')
plt.title('Distribution of Gene Expression')
plt.show()
```

```python
# Line plot (e.g., expression across samples)
import matplotlib.pyplot as plt
samples = ['S1', 'S2', 'S3', 'S4']
expr = [10, 15, 12, 18]
plt.plot(samples, expr, marker='o')
plt.xlabel('Sample')
plt.ylabel('Expression')
plt.title('Expression Across Samples')
plt.legend(['Gene X'])
plt.show()
```

### R: base graphics and ggplot2

```r
# Scatter plot (base R)
plot(df$expression, df$pvalue,
     xlab = "Expression (TPM)",
     ylab = "P-value",
     main = "Gene Expression vs Significance")
```

```r
# Histogram
hist(df$expression, breaks = 20,
     xlab = "Expression",
     main = "Distribution of Gene Expression")
```

```r
# With ggplot2 (if installed):
# library(ggplot2)
# ggplot(df, aes(x=expression, y=pvalue)) + geom_point() + labs(title="Expression vs P-value")
```

**Practice:** Create a scatter plot with labeled axes and a title.

---

## 6. Writing Functions

### Python

```python
def calculate_mean(numbers):
    """Return the average of a list of numbers."""
    return sum(numbers) / len(numbers)

expression_values = [10.2, 5.1, 0.3, 22.0, 8.7]
avg = calculate_mean(expression_values)
print(avg)
```

```python
def calculate_median(numbers):
    """Return the median (middle value when sorted)."""
    sorted_nums = sorted(numbers)
    n = len(sorted_nums)
    mid = n // 2
    if n % 2 == 1:
        return sorted_nums[mid]
    return (sorted_nums[mid - 1] + sorted_nums[mid]) / 2

print(calculate_median([1, 3, 5, 7, 9]))  # 5
```

### R

```r
calculate_mean <- function(numbers) {
  return(sum(numbers) / length(numbers))
}

expression_values <- c(10.2, 5.1, 0.3, 22.0, 8.7)
avg <- calculate_mean(expression_values)
print(avg)
```

```r
calculate_median <- function(numbers) {
  sorted_nums <- sort(numbers)
  n <- length(sorted_nums)
  mid <- ceiling(n / 2)
  if (n %% 2 == 1) {
    return(sorted_nums[mid])
  }
  return((sorted_nums[mid - 1] + sorted_nums[mid]) / 2)
}

print(calculate_median(c(1, 3, 5, 7, 9)))  # 5
```

**Practice:** Write a `calculate_mean` function and test it with `[2, 4, 6, 8, 10]`.

---

## 7. Conditionals (if / elif / else)

### Python

```python
# Syntax: indentation matters!
age = 25
if age >= 18:
    print("Adult")
elif age >= 13:
    print("Teen")
else:
    print("Child")
```

```python
# Comparison operators: ==, !=, <, >, <=, >=
pvalue = 0.03
if pvalue < 0.05:
    print("Significant")
else:
    print("Not significant")
```

```python
# Nested conditionals
fc = 2.5
pval = 0.01
if pval < 0.05:
    if fc > 2:
        print("Significant and upregulated")
    else:
        print("Significant but not strongly upregulated")
```

### R

```r
age <- 25
if (age >= 18) {
  print("Adult")
} else if (age >= 13) {
  print("Teen")
} else {
  print("Child")
}
```

```r
# R uses double-equals for equality
pvalue <- 0.03
if (pvalue < 0.05) {
  print("Significant")
} else {
  print("Not significant")
}
```

**Practice:** Write code that prints "High" if expression > 10, "Low" if < 1, else "Medium".

---

## 8. Loops (for, while)

### Python

```python
# For loop: squares of 1 to 5
for i in range(1, 6):
    print(i ** 2)
```

```python
# Accumulate: sum expression values
total = 0
for val in [10, 15, 20]:
    total += val
print(total)
```

```python
# While loop: count until threshold
count = 0
while count < 5:
    print(count)
    count += 1
```

### R

```r
# For loop
for (i in 1:5) {
  print(i^2)
}
```

```r
# Vectorized alternative (R idiom)
print((1:5)^2)
```

```r
# While loop
count <- 0
while (count < 5) {
  print(count)
  count <- count + 1
}
```

**Practice:** Use a loop to print the cube of each number from 1 to 4.

---

## 9. Dictionaries (Python) / Named Lists (R)

### Python: Dictionaries

Dictionaries map keys to values. Perfect for gene → expression lookups.

```python
genes = {'TP53': 100, 'BRCA1': 50, 'MYC': 75}
print(genes['BRCA1'])       # 50
print(genes.get('BRCA1'))   # 50 (safer)
print(genes.get('UNKNOWN', 0))  # 0 (default if missing)
```

```python
# Iterate over keys and values
for gene, expr in genes.items():
    print(f"{gene}: {expr} TPM")
```

```python
# Modify
genes['TP53'] = 120
genes['EGFR'] = 90
```

### R: Named Vectors / Lists

```r
# Named list (or named vector)
genes <- list(TP53 = 100, BRCA1 = 50, MYC = 75)
print(genes$BRCA1)          # 50
print(genes[["BRCA1"]])     # Alternative syntax
```

```r
# Named vector
genes_vec <- c(TP53 = 100, BRCA1 = 50, MYC = 75)
print(genes_vec["BRCA1"])
```

```r
# Iterate
for (gene in names(genes)) {
  print(paste(gene, ":", genes[[gene]], "TPM"))
}
```

**Practice:** Create a dictionary/list with 3 gene names and expression values. Retrieve one value and add a new gene.

---

## 10. Error Types and Debugging

### Python: Common Errors

| Error       | Cause                                  |
|-------------|----------------------------------------|
| IndexError  | Accessing `my_list[10]` when len = 5   |
| KeyError    | Accessing missing dict key             |
| TypeError   | Wrong type (e.g., int + string)        |
| NameError   | Variable not defined                   |

```python
# IndexError
my_list = [1, 2, 3, 4, 5]
# my_list[10]  # IndexError: list index out of range

# KeyError
# genes['UNKNOWN']  # KeyError if key doesn't exist
genes.get('UNKNOWN', 0)  # Safe: returns 0

# TypeError
# "hello" + 5  # TypeError: can only concatenate str + str

# NameError
# print(undefined_var)  # NameError: name 'undefined_var' is not defined
```

**Reading error messages:** The last line shows the error type. The traceback shows where it occurred.

### R: Common Errors and Packages

R has similar concepts. For subscript out of bounds, you get `subscript out of bounds`. For packages:

```r
# Install packages from CRAN
install.packages("dplyr")   # One-time install
library(dplyr)              # Load for session

# Common R errors:
# - "subscript out of bounds" (like IndexError)
# - "object 'x' not found" (like NameError)
# - "could not find function" (package not loaded)
```

**Debugging strategies:**
1. Read the error message carefully
2. Add `print()` statements to inspect values
3. Check indices and lengths
4. Use a debugger or IDE breakpoints

---

## Practice Exercises

### Exercise 1: Combine Syntax and Data Structures

**Task:** Create a list/vector of 5 gene expression values. Compute the mean using a loop (without built-in mean). Print the result.

<details>
<summary>Python solution</summary>

```python
expr = [10.2, 5.1, 0.3, 22.0, 8.7]
total = 0
for val in expr:
    total += val
mean = total / len(expr)
print(mean)
```
</details>

<details>
<summary>R solution</summary>

```r
expr <- c(10.2, 5.1, 0.3, 22.0, 8.7)
total <- 0
for (val in expr) {
  total <- total + val
}
mean_val <- total / length(expr)
print(mean_val)
```
</details>

---

### Exercise 2: Data Manipulation and Conditionals

**Task:** Load a CSV with columns `gene`, `expression`, `pvalue`. Filter rows where `pvalue < 0.05` and `expression > 5`. Print the number of such genes.

<details>
<summary>Python solution</summary>

```python
import pandas as pd
df = pd.read_csv("your_data.csv")
filtered = df[(df['pvalue'] < 0.05) & (df['expression'] > 5)]
print(len(filtered))
```
</details>

<details>
<summary>R solution</summary>

```r
df <- read.csv("your_data.csv")
filtered <- df[df$pvalue < 0.05 & df$expression > 5, ]
print(nrow(filtered))
```
</details>

---

### Exercise 3: Function and List/Vector

**Task:** Write `calculate_median` and test with `[3, 1, 4, 1, 5]`. Expected output: 3.

<details>
<summary>Python solution</summary>

```python
def calculate_median(numbers):
    sorted_nums = sorted(numbers)
    n = len(sorted_nums)
    mid = n // 2
    if n % 2 == 1:
        return sorted_nums[mid]
    return (sorted_nums[mid - 1] + sorted_nums[mid]) / 2
print(calculate_median([3, 1, 4, 1, 5]))
```
</details>

<details>
<summary>R solution</summary>

```r
calculate_median <- function(numbers) {
  sorted_nums <- sort(numbers)
  n <- length(sorted_nums)
  mid <- ceiling(n / 2)
  if (n %% 2 == 1) {
    return(sorted_nums[mid])
  }
  (sorted_nums[mid - 1] + sorted_nums[mid]) / 2
}
print(calculate_median(c(3, 1, 4, 1, 5)))
```
</details>

---

### Exercise 4: Dictionary/List and Loop

**Task:** Given `genes = {'A': 10, 'B': 20, 'C': 5}`, print each gene name and whether its expression is "high" (>15) or "low" (≤15).

<details>
<summary>Python solution</summary>

```python
genes = {'A': 10, 'B': 20, 'C': 5}
for gene, expr in genes.items():
    if expr > 15:
        print(f"{gene}: high")
    else:
        print(f"{gene}: low")
```
</details>

<details>
<summary>R solution</summary>

```r
genes <- list(A = 10, B = 20, C = 5)
for (gene in names(genes)) {
  expr <- genes[[gene]]
  if (expr > 15) {
    print(paste(gene, ": high"))
  } else {
    print(paste(gene, ": low"))
  }
}
```
</details>

---

### Exercise 5: Full Workflow

**Task:** Create a small in-memory dataset (3 genes, expression and pvalue). Filter significant genes (p < 0.05), sort by expression descending, and plot expression vs pvalue.

<details>
<summary>Python solution</summary>

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame({
    'gene': ['TP53', 'BRCA1', 'MYC'],
    'expression': [100, 50, 75],
    'pvalue': [0.01, 0.04, 0.10]
})
filtered = df[df['pvalue'] < 0.05].sort_values('expression', ascending=False)
plt.scatter(df['expression'], df['pvalue'])
plt.xlabel('Expression')
plt.ylabel('P-value')
plt.title('Gene Expression vs P-value')
plt.show()
```
</details>

<details>
<summary>R solution</summary>

```r
df <- data.frame(
  gene = c("TP53", "BRCA1", "MYC"),
  expression = c(100, 50, 75),
  pvalue = c(0.01, 0.04, 0.10)
)
filtered <- df[df$pvalue < 0.05, ]
filtered <- filtered[order(filtered$expression, decreasing = TRUE), ]
plot(df$expression, df$pvalue, xlab = "Expression", ylab = "P-value",
     main = "Gene Expression vs P-value")
```
</details>

---

## Additional Resources

- **Python:** [Learn Python](https://www.learnpython.org/)
- **R:** [Swirl Stats](https://swirlstats.com/students.html)
- **Both:** [DataCamp](https://www.datacamp.com/) (free tier)
- **Python for Biology:** [Rosalind](https://rosalind.info/) (bioinformatics problems)
- **R for Data Science:** [R4DS](https://r4ds.hadley.nz/) (free online book)
