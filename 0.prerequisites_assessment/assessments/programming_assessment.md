# Programming Assessment

**Time Limit:** 30 minutes  
**Passing Score:** 7/10 correct  
**Language:** Choose R or Python (whichever you'll use for the course)

---

## Instructions

1. Answer all questions for your chosen language
2. For coding questions, write the code (or pseudocode if unsure)
3. Check your answers against the answer key (bottom of file)
4. If you score < 70%, review the programming primer before starting Course 1

---

## Python Assessment

### Question 1 (Basic Syntax)
What will this code print?
```python
x = 5
y = 10
print(x + y * 2)
```
**A)** 30  
**B)** 25  
**C)** 20  
**D)** Error  

**Your answer:** _____

---

### Question 2 (Lists)
Write code to create a list of numbers from 1 to 10 and print only the even numbers.

**Your code:**
```python
# Write your answer here
```

---

### Question 3 (Reading Data)
Which library is most commonly used to read CSV files in Python for data analysis?

**A)** csv  
**B)** pandas  
**C)** numpy  
**D)** matplotlib  

**Your answer:** _____

---

### Question 4 (Data Manipulation)
Given a pandas DataFrame `df` with a column 'score', write code to filter rows where score > 50.

**Your code:**
```python
# Write your answer here
```

---

### Question 5 (Plotting)
What library would you use to create a basic scatter plot in Python?

**A)** pandas  
**B)** numpy  
**C)** matplotlib or seaborn  
**D)** scipy  

**Your answer:** _____

---

### Question 6 (Functions)
Write a function called `calculate_mean` that takes a list of numbers and returns their average.

**Your code:**
```python
# Write your answer here
```

---

### Question 7 (Conditionals)
What will this code print?
```python
age = 25
if age >= 18:
    print("Adult")
elif age >= 13:
    print("Teen")
else:
    print("Child")
```

**Your answer:** _____

---

### Question 8 (Loops)
Write code to print the square of each number from 1 to 5.

**Your code:**
```python
# Write your answer here
```

---

### Question 9 (Dictionaries)
Given a dictionary `genes = {'TP53': 100, 'BRCA1': 50, 'MYC': 75}`, write code to get the expression value of 'BRCA1'.

**Your code:**
```python
# Write your answer here
```

---

### Question 10 (Error Handling)
What type of error would occur if you try to access `my_list[10]` when `my_list` only has 5 elements?

**A)** TypeError  
**B)** IndexError  
**C)** ValueError  
**D)** KeyError  

**Your answer:** _____

---

## R Assessment

### Question 1 (Basic Syntax)
What will this code print?
```r
x <- 5
y <- 10
print(x + y * 2)
```
**A)** 30  
**B)** 25  
**C)** 20  
**D)** Error  

**Your answer:** _____

---

### Question 2 (Vectors)
Write code to create a vector of numbers from 1 to 10 and print only the even numbers.

**Your code:**
```r
# Write your answer here
```

---

### Question 3 (Reading Data)
Which function is commonly used to read CSV files in R?

**A)** read.csv()  
**B)** load.csv()  
**C)** import.csv()  
**D)** get.csv()  

**Your answer:** _____

---

### Question 4 (Data Manipulation)
Given a data frame `df` with a column 'score', write code to filter rows where score > 50.

**Your code:**
```r
# Write your answer here
```

---

### Question 5 (Plotting)
What is the basic plotting system in R called?

**A)** matplotlib  
**B)** base R graphics or ggplot2  
**C)** seaborn  
**D)** plotly  

**Your answer:** _____

---

### Question 6 (Functions)
Write a function called `calculate_mean` that takes a vector of numbers and returns their average.

**Your code:**
```r
# Write your answer here
```

---

### Question 7 (Conditionals)
What will this code print?
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

**Your answer:** _____

---

### Question 8 (Loops)
Write code to print the square of each number from 1 to 5.

**Your code:**
```r
# Write your answer here
```

---

### Question 9 (Lists)
Given a list `genes <- list(TP53=100, BRCA1=50, MYC=75)`, write code to get the expression value of 'BRCA1'.

**Your code:**
```r
# Write your answer here
```

---

### Question 10 (Packages)
What function do you use to install packages from CRAN in R?

**A)** install()  
**B)** install.packages()  
**C)** load.package()  
**D)** get.package()  

**Your answer:** _____

---

## Scoring

Count your correct answers: _____ / 10

**Interpretation:**
- **9-10:** Excellent! You're well-prepared.
- **7-8:** Good foundation. Review any missed topics.
- **5-6:** Borderline. Study the programming primer, then retake.
- **0-4:** Need significant study. Complete programming primer before retaking.

---

## Answer Key

### Python Answers

1. **B) 25** (Order of operations: 10 * 2 = 20, then 5 + 20 = 25)

2. **Correct solutions include:**
```python
numbers = list(range(1, 11))
for num in numbers:
    if num % 2 == 0:
        print(num)
# OR
numbers = list(range(1, 11))
evens = [n for n in numbers if n % 2 == 0]
print(evens)
```

3. **B) pandas**

4. **Correct solutions include:**
```python
filtered_df = df[df['score'] > 50]
# OR
filtered_df = df.query('score > 50')
```

5. **C) matplotlib or seaborn**

6. **Correct solutions include:**
```python
def calculate_mean(numbers):
    return sum(numbers) / len(numbers)
# OR
import statistics
def calculate_mean(numbers):
    return statistics.mean(numbers)
# OR
import numpy as np
def calculate_mean(numbers):
    return np.mean(numbers)
```

7. **"Adult"**

8. **Correct solutions include:**
```python
for i in range(1, 6):
    print(i ** 2)
# OR
for i in range(1, 6):
    print(i * i)
```

9. **Correct solutions include:**
```python
expression = genes['BRCA1']
# OR
expression = genes.get('BRCA1')
```

10. **B) IndexError**

---

### R Answers

1. **B) 25** (Order of operations: 10 * 2 = 20, then 5 + 20 = 25)

2. **Correct solutions include:**
```r
numbers <- 1:10
evens <- numbers[numbers %% 2 == 0]
print(evens)
# OR
for (num in 1:10) {
  if (num %% 2 == 0) {
    print(num)
  }
}
```

3. **A) read.csv()**

4. **Correct solutions include:**
```r
filtered_df <- df[df$score > 50, ]
# OR (using dplyr)
library(dplyr)
filtered_df <- filter(df, score > 50)
# OR
filtered_df <- subset(df, score > 50)
```

5. **B) base R graphics or ggplot2**

6. **Correct solutions include:**
```r
calculate_mean <- function(numbers) {
  return(mean(numbers))
}
# OR
calculate_mean <- function(numbers) {
  sum(numbers) / length(numbers)
}
```

7. **"Adult"**

8. **Correct solutions include:**
```r
for (i in 1:5) {
  print(i^2)
}
# OR
for (i in 1:5) {
  print(i * i)
}
```

9. **Correct solutions include:**
```r
expression <- genes$BRCA1
# OR
expression <- genes[["BRCA1"]]
# OR
expression <- genes[[2]]
```

10. **B) install.packages()**

---

## Next Steps

### If You Passed (≥ 7/10):
✅ Move to next assessment: `command_line_assessment.md`

### If You Need More Practice:
📚 Review: `../resources/programming_primer.md`
- Focus on areas where you missed questions
- Complete practice exercises
- Retake this assessment when ready

### Additional Practice Resources:
- **Python:** https://www.learnpython.org/
- **R:** https://swirlstats.com/students.html
- **Both:** https://www.datacamp.com/ (free tier)
