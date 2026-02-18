# Command Line Assessment

**Time Limit:** 20 minutes  
**Passing Score:** 7/10 correct  

---

## Instructions

1. Answer all questions using standard Unix/Linux/macOS commands
2. For command questions, write the exact command(s) you would type
3. Check your answers against the answer key (bottom of file)
4. If you score < 70%, review the command line primer before starting Course 1

---

### Question 1 (Navigation)
What does the `pwd` command do?

**A)** Print working directory  
**B)** Create a new directory  
**C)** List files in current directory  
**D)** Change to parent directory  

**Your answer:**

---

### Question 2 (File Operations)
Write the command to create a directory called `results` and copy the file `data.txt` into it.

**Your command:**
```bash
# Write your answer here

```

---

### Question 3 (Viewing File Contents)
Which command would you use to view only the first 15 lines of a large file?

**A)** `tail`  
**B)** `cat`  
**C)** `head`  
**D)** `less`  

**Your answer:**

---

### Question 4 (Pipes and Redirects)
Write the command to run `ls -l` and save its output to a file named `file_list.txt` (overwriting if it exists).

**Your command:**
```bash
# Write your answer here

```

---

### Question 5 (File Permissions)
In `ls -l` output, what does `drwxr-xr-x` mean for the first character `d`?

**A)** The file is executable  
**B)** The entry is a directory  
**C)** The file is a symbolic link  
**D)** The file has special permissions  

**Your answer:**

---

### Question 6 (grep and Pattern Matching)
Write the command to search for lines containing the word "error" (case-insensitive) in the file `server.log`.

**Your command:**
```bash
# Write your answer here

```

---

### Question 7 (Environment Variables and PATH)
How do you display the value of the `PATH` environment variable?

**A)** `show PATH`  
**B)** `echo $PATH`  
**C)** `print PATH`  
**D)** `PATH`  

**Your answer:**

---

### Question 8 (Process Management)
Write the command to terminate a process with PID 12345.

**Your command:**
```bash
# Write your answer here

```

---

### Question 9 (Compression and Archives)
What does the command `tar -xzf archive.tar.gz` do?

**A)** Create a compressed archive  
**B)** Extract the contents of a gzipped tar archive  
**C)** List the contents without extracting  
**D)** Compress a single file  

**Your answer:**

---

### Question 10 (SSH and Remote Connections)
Write the command to connect via SSH to a remote server as user `scientist` on host `lab.university.edu`.

**Your command:**
```bash
# Write your answer here

```

---

## Scoring

Count your correct answers: _____ / 10

**Interpretation:**
- **9-10:** Excellent! You're well-prepared for the command line.
- **7-8:** Good foundation. Review any missed topics.
- **5-6:** Borderline. Study the command line primer, then retake.
- **0-4:** Need significant study. Complete command line primer before retaking.

---

## Answer Key

1. **A) Print working directory** — `pwd` displays the full path of your current working directory.

2. **Correct solutions include:**
```bash
mkdir results
cp data.txt results/
# OR in one line:
mkdir results && cp data.txt results/
```

3. **C) head** — `head` displays the first lines of a file; use `head -n 15 filename` for the first 15 lines. (`less` could work interactively but `head` directly shows the first N lines.)

4. **Correct solutions include:**
```bash
ls -l > file_list.txt
```

5. **B) The entry is a directory** — The first character in the permissions string indicates file type: `d` = directory, `-` = regular file, `l` = symbolic link.

6. **Correct solutions include:**
```bash
grep -i error server.log
# OR
grep -i "error" server.log
```

7. **B) echo $PATH** — In bash/zsh, you use `echo $VARIABLE` to display environment variable values. The `$` dereferences the variable.

8. **Correct solutions include:**
```bash
kill 12345
# OR (for forceful termination)
kill -9 12345
```

9. **B) Extract the contents of a gzipped tar archive** — `tar -xzf` means: `-x` extract, `-z` gzip, `-f` file. It decompresses and extracts the archive.

10. **Correct solutions include:**
```bash
ssh scientist@lab.university.edu
```

---

## Next Steps

### If You Passed (≥ 7/10):
✅ Move to next assessment: `biology_stats_assessment.md`

### If You Need More Practice:
📚 Review: `../resources/command_line_primer.md`
- Focus on areas where you missed questions
- Practice commands in a terminal
- Retake this assessment when ready
