# Git Assessment

**Time Limit:** 20 minutes  
**Passing Score:** 6/10 correct (60%)

---

## Instructions

1. Answer all questions
2. For command questions, write the exact Git command(s)
3. Check your answers against the answer key (bottom of file)
4. If you score < 60%, review the Git primer before starting Course 1

---

### Question 1 (Version Control Basics)
What is Git primarily used for?

**A)** A code editor  
**B)** Version control—tracking changes to files over time  
**C)** A programming language  
**D)** A web browser  

**Your answer:**

---

### Question 2 (Clone)
Write the command to clone a repository from the URL `https://github.com/user/repo.git`.

**Your command:**
```bash
# Write your answer here
```

---

### Question 3 (Status)
What does `git status` show?

**A)** The commit history  
**B)** The current branch and which files are modified, staged, or untracked  
**C)** A list of remote repositories  
**D)** The diff between your code and the remote  

**Your answer:**

---

### Question 4 (Staging)
Write the command to stage all changed files for the next commit.

**Your command:**
```bash
# Write your answer here
```

---

### Question 5 (Commit)
Write the command to commit your staged changes with the message "Fix typo in README".

**Your command:**
```bash
# Write your answer here
```

---

### Question 6 (Push and Pull)
What is the difference between `git push` and `git pull`?

**A)** Both do the same thing  
**B)** `git push` sends your commits to the remote; `git pull` fetches and merges changes from the remote  
**C)** `git pull` sends your commits; `git push` fetches from the remote  
**D)** They are used only for creating branches  

**Your answer:**

---

### Question 7 (Branches)
What is a branch in Git?

**A)** A backup copy of the entire repository  
**B)** An independent line of development that allows you to work on features or fixes without affecting the main codebase  
**C)** A remote server that stores the repository  
**D)** A way to delete old commits  

**Your answer:**

---

### Question 8 (Create and Switch Branches)
Write the command to create a new branch called `feature-login` and switch to it.

**Your command:**
```bash
# Write your answer here
```

---

### Question 9 (.gitignore)
What is the purpose of a `.gitignore` file?

**A)** To ignore all files in the repository  
**B)** To specify which files or patterns Git should not track (e.g., build outputs, secrets)  
**C)** To delete files from the repository  
**D)** To lock the repository from further changes  

**Your answer:**

---

### Question 10 (Merge Conflicts)
When do merge conflicts occur?

**A)** Every time you run `git merge`  
**B)** When Git cannot automatically reconcile different changes to the same lines of the same file  
**C)** Only when working on different branches  
**D)** When the remote repository is offline  

**Your answer:**

---

## Scoring

Count your correct answers: _____ / 10

**Interpretation:**
- **9-10:** Excellent! You're well-prepared.
- **7-8:** Good foundation. Review any missed topics.
- **5-6:** Borderline. Study the Git primer, then retake.
- **0-4:** Need significant study. Complete the Git primer before retaking.

---

## Answer Key

1. **B) Version control—tracking changes to files over time** — Git is a distributed version control system that tracks changes, preserves history, and enables collaboration.

2. **Correct command:**
```bash
git clone https://github.com/user/repo.git
```

3. **B) The current branch and which files are modified, staged, or untracked** — `git status` shows your working tree state: current branch, modified files, staged changes, and untracked files.

4. **Correct command:**
```bash
git add .
# OR to add specific files:
git add file1.txt file2.txt
```

5. **Correct command:**
```bash
git commit -m "Fix typo in README"
```

6. **B) `git push` sends your commits to the remote; `git pull` fetches and merges changes from the remote** — Push uploads your local commits; pull downloads and merges remote changes into your branch.

7. **B) An independent line of development that allows you to work on features or fixes without affecting the main codebase** — Branches let you isolate work and merge it when ready.

8. **Correct commands include:**
```bash
git checkout -b feature-login
# OR (Git 2.23+):
git switch -c feature-login
```

9. **B) To specify which files or patterns Git should not track (e.g., build outputs, secrets)** — `.gitignore` prevents sensitive or generated files from being committed.

10. **B) When Git cannot automatically reconcile different changes to the same lines of the same file** — Conflicts arise when two branches modify the same lines; you must manually resolve them.

---

## Next Steps

### If You Passed (≥ 6/10):
✅ Move to next assessment: `environment_test.md`

### If You Need More Practice:
📚 Review: `../resources/git_primer.md`
- Focus on areas where you missed questions
- Complete practice exercises
- Retake this assessment when ready
