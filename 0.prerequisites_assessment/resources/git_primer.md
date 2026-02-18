# Git Primer for Bioinformatics

**Estimated Time:** 3-5 hours (self-paced)  
**Goal:** Build Git skills to pass the Git Assessment (6/10)

---

## How to Use This Primer

- Install Git first if you haven't
- Follow along by typing every command
- Create a practice repo as you go

---

## Topics

### 1. What Is Version Control?

**Why version control matters**

When you analyze data, write scripts, or develop pipelines, your work evolves over time. Without version control, you end up with files like `analysis_v1.R`, `analysis_v2.R`, and `final_v2_FINAL_really_final.R`—and no clear record of what changed or why.

**The problem: final_v2_FINAL_really_final.R**

This chaos leads to:
- Lost work when you overwrite the wrong file
- Confusion about which version is current
- Inability to undo mistakes
- Nightmares when collaborating

**Git as distributed version control**

Git tracks every change you make, lets you revert to any previous state, and supports collaboration. Unlike centralized systems, every copy of a Git repository contains the full history—you can work offline and sync later.

**Repository concept**

A **repository** (or "repo") is a folder where Git tracks all changes. It contains your project files plus a hidden `.git` directory that stores the complete history.

---

### 2. Setting Up Git

**Installing Git**

- **Linux:** `sudo apt install git` (Debian/Ubuntu) or `sudo dnf install git` (Fedora)
- **macOS:** `brew install git` or install Xcode Command Line Tools
- **Windows:** Download from [git-scm.com](https://git-scm.com)

**git config (name, email)**

Before your first commit, configure your identity:

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

Use the same email as your GitHub/GitLab account for commits to link correctly.

**Checking your setup**

```bash
git --version
git config --list
```

---

### 3. Creating and Cloning Repositories

**git init (new repo)**

To start version control in an existing folder:

```bash
cd my_project
git init
```

**git clone (existing repo)**

To get a copy of an existing repository (e.g., from GitHub):

```bash
git clone https://github.com/username/repo-name.git
cd repo-name
```

**The .git directory**

The `.git` folder contains all version history. Never delete it manually unless you want to destroy the repo. It is hidden by default.

---

### 4. The Git Workflow: Edit → Stage → Commit

**Working directory → staging area → repository**

1. **Working directory:** Your actual files. Edit them as usual.
2. **Staging area (index):** A "holding zone" for changes you want to include in the next commit.
3. **Repository:** The committed history. Once committed, changes are permanently recorded.

```
[Working Dir] --git add--> [Staging Area] --git commit--> [Repository]
```

**Diagram/explanation of the three states**

- **Modified:** You changed a file but haven't staged it.
- **Staged:** You added the file; it's ready to commit.
- **Committed:** The change is safely stored in the repository.

**git status (checking state)**

```bash
git status
```

Run this often. It tells you which files are modified, staged, or untracked.

---

### 5. Staging Changes

**git add (specific files, all files)**

```bash
git add script.R              # Stage one file
git add script.R utils.R      # Stage multiple files
git add src/                  # Stage entire directory
```

**git add . vs git add -A**

- `git add .` — Stages new and modified files in the current directory and subdirectories (but behavior can vary with Git version).
- `git add -A` — Stages all changes everywhere: new, modified, and deleted files.

For clarity, many prefer `git add -A` or `git add --all`.

**Unstaging: git restore --staged**

To remove a file from the staging area without discarding changes:

```bash
git restore --staged filename.R
```

(Older Git: `git reset HEAD filename.R`)

---

### 6. Committing Changes

**git commit -m "message"**

```bash
git commit -m "Add QC step to RNA-seq pipeline"
```

**Writing good commit messages**

- Use the imperative: "Add feature" not "Added feature"
- Be specific: "Fix sample ID parsing in metadata" not "Fix bug"
- Keep the first line under ~50 characters; add detail in the body if needed

**Atomic commits (one logical change per commit)**

Each commit should represent one logical change. Don't mix unrelated edits (e.g., fixing a typo and adding a new function) in one commit. This makes history easier to understand and revert.

**git log to view history**

```bash
git log
git log --oneline          # Compact one-line view
git log -n 5               # Last 5 commits
```

---

### 7. Working with Remotes

**What is a remote (origin)**

A **remote** is a version of your repository hosted elsewhere (e.g., GitHub, GitLab). `origin` is the default name for the repository you cloned from or first pushed to.

**git push (upload changes)**

```bash
git push origin main
```

Uploads your local commits to the remote. Replace `main` with your branch name if different (e.g., `master`).

**git pull (download changes)**

```bash
git pull origin main
```

Downloads remote changes and merges them into your current branch. Run this before you start working and before pushing.

**git fetch vs git pull**

- `git fetch` — Downloads remote changes but does not merge. Lets you inspect before integrating.
- `git pull` — Fetches and merges in one step. Convenient but can surprise you with merge commits.

For beginners, `git pull` before you start work is usually enough.

---

### 8. Branches

**What branches are and why they matter**

A **branch** is a parallel line of development. You can work on a new feature without affecting the main code until you merge. Essential for collaboration and experimentation.

**git branch (list/create)**

```bash
git branch                    # List branches
git branch feature-qc         # Create new branch
```

**git checkout / git switch**

```bash
git checkout feature-qc       # Switch to branch (older syntax)
git switch feature-qc         # Switch to branch (newer, clearer)
```

**Working on feature branches**

1. Create and switch: `git switch -c feature-name`
2. Make changes and commit
3. Switch back to main: `git switch main`
4. Merge: `git merge feature-name`

---

### 9. Merging and Conflicts

**git merge basics**

```bash
git switch main
git merge feature-branch
```

Brings changes from `feature-branch` into `main`.

**What merge conflicts look like**

When Git cannot automatically combine changes, you get a **merge conflict**. Affected files contain markers:

```
<<<<<<< HEAD
current version on your branch
=======
incoming version from the other branch
>>>>>>> feature-branch
```

**How to resolve: the <<<< ==== >>>> markers**

1. Open the file in an editor.
2. Find the `<<<<<<<`, `=======`, `>>>>>>>` blocks.
3. Decide what the final content should be (keep one version, combine both, or write something new).
4. Delete the conflict markers.
5. Stage and commit: `git add file.R` then `git commit`.

**Best practices to avoid conflicts**

- Pull frequently so your branch stays close to main.
- Communicate with collaborators about who edits what.
- Keep commits small and focused.
- Use `.gitignore` so generated files don't clutter the repo.

---

### 10. .gitignore

**Purpose and syntax**

A `.gitignore` file tells Git which files or folders to ignore. One pattern per line. `#` starts a comment.

**Common patterns for bioinformatics**

```
# Large data files
*.fastq
*.fastq.gz
*.bam
*.bai
*.vcf
*.vcf.gz

# Environment and secrets
.env
.env.local
*.key

# Output directories
output/
results/
plots/
*.RData

# OS files
.DS_Store
Thumbs.db

# IDE
.idea/
.vscode/
*.swp
```

**Template .gitignore example**

Create a `.gitignore` in your repo root:

```
# Data
*.fastq
*.fastq.gz
*.bam
*.bai

# R
.RData
.Rhistory

# Python
__pycache__/
*.pyc
venv/
.venv/

# Outputs
output/
results/
```

---

### 11. Essential Commands Quick Reference

| Command | Description |
|---------|-------------|
| `git init` | Initialize a new repository |
| `git clone <url>` | Clone an existing repository |
| `git status` | Check working tree status |
| `git add <file>` | Stage changes |
| `git commit -m "msg"` | Commit staged changes |
| `git log` | View commit history |
| `git pull` | Fetch and merge from remote |
| `git push` | Push commits to remote |
| `git branch` | List branches |
| `git switch <branch>` | Switch branches |
| `git merge <branch>` | Merge a branch into current |
| `git restore --staged <file>` | Unstage a file |

---

### 12. Best Practices for Bioinformatics

- **Commit analysis scripts, not large data** — Use `.gitignore` for raw data, alignment files, and large outputs. Store scripts, configs, and small metadata.
- **Use meaningful commit messages** — Future you (and collaborators) will thank you.
- **Keep repos organized** — Use clear folder structure (`scripts/`, `data/`, `output/`, `docs/`).
- **Tag important versions** — Use `git tag v1.0` for manuscript-ready or release versions. Easy to find later: `git checkout v1.0`.

---

## Practice Exercises

### Exercise 1: Create a repository

1. Create a folder: `mkdir git_practice && cd git_practice`
2. Initialize: `git init`
3. Create a file: `echo "# My Analysis" > README.md`
4. Stage and commit: `git add README.md` then `git commit -m "Initial commit"`

### Exercise 2: Make multiple commits

1. Edit `README.md` (add a line)
2. Run `git status`
3. Stage and commit with a descriptive message
4. Run `git log --oneline` to see your history

### Exercise 3: Create a branch

1. Create and switch: `git switch -c add-scripts`
2. Create `scripts/hello.R` with a simple R script
3. Stage and commit
4. Switch back: `git switch main`
5. Verify the file is not on main: `ls scripts/` (or it won't exist yet)

### Exercise 4: Merge a branch

1. From `main`, run: `git merge add-scripts`
2. Check that `scripts/hello.R` now exists on main
3. Run `git log` to see the merge

### Exercise 5: Resolve a merge conflict

1. On `main`, edit `README.md` (e.g., change line 1)
2. Commit the change
3. Switch to `add-scripts` and edit the same line of `README.md` differently
4. Commit
5. Switch to `main` and run `git merge add-scripts`
6. Open `README.md`, resolve the conflict markers, then `git add README.md` and `git commit`

---

## Additional Resources

- [GitHub Learning Lab](https://lab.github.com/) — Free interactive Git courses
- [Git Immersion](https://gitimmersion.com/) — Step-by-step Git tutorial
- [Pro Git Book](https://git-scm.com/book/en/v2) — Free, comprehensive reference
- [Git Cheat Sheet](https://education.github.com/git-cheat-sheet-education.pdf) — Quick reference PDF
