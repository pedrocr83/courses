# Command Line Primer for Bioinformatics

**Estimated Time:** 5-10 hours (self-paced)  
**Goal:** Build command line skills to pass the Command Line Assessment (7/10)

---

## How to Use This Primer

- Work through each section in order
- Type every command yourself (don't just read)
- Complete practice exercises

---

## 1. Getting Started with the Terminal

### What is a Shell (bash)?

The **shell** is a program that interprets your typed commands and talks to the operating system. **Bash** (Bourne Again Shell) is the default on most Linux systems and macOS. It's your interface to run programs, navigate files, and automate tasks.

### Opening a Terminal

- **Linux:** `Ctrl+Alt+T` or search for "Terminal"
- **macOS:** `Cmd+Space` → type "Terminal" → Enter
- **Windows:** Use Windows Subsystem for Linux (WSL) or Git Bash

### The Prompt

When you open a terminal, you see a **prompt**—text waiting for your input. It often shows your username, hostname, and current directory:

```
user@machine:~$
```

Everything after the `$` (or `%` on some systems) is where you type commands.

---

## 2. Navigation (cd, ls, pwd)

### Directory Structure (Tree Metaphor)

The filesystem is a tree of **directories** (folders) and **files**. The root `/` is the top. Your home directory is usually `/home/username` (Linux) or `/Users/username` (macOS).

### Key Commands

| Command | Purpose |
|---------|---------|
| `pwd` | **P**rint **w**orking **d**irectory—shows where you are |
| `ls` | List files and directories in the current location |
| `cd` | Change directory |

### Absolute vs Relative Paths

- **Absolute:** Starts from root (`/home/user/projects`)
- **Relative:** Starts from current directory (`../data` or `results/`)

### Special Paths

| Symbol | Meaning |
|--------|---------|
| `~` | Your home directory |
| `.` | Current directory |
| `..` | Parent directory |

```bash
cd ~              # Go to home directory
cd ..             # Go up one level
cd ../data        # Go to sibling 'data' directory
```

### ls with Useful Flags

```bash
ls -l             # Long format (permissions, owner, size, date)
ls -a             # Show hidden files (those starting with .)
ls -lh            # Human-readable sizes (K, M, G)
ls -la            # Combine: long + hidden
```

---

## 3. File and Directory Operations

### Creating and Removing Directories

```bash
mkdir results              # Create directory 'results'
mkdir -p analysis/output   # Create nested dirs; -p creates parents as needed
rmdir empty_dir            # Remove empty directory only
```

### File Operations

```bash
touch data.txt             # Create empty file or update timestamp
cp data.txt backup/        # Copy file to backup/
cp -r project/ archive/    # Copy directory recursively (-r)
mv data.txt renamed.txt    # Move (rename) file
mv old/ new/               # Rename directory
rm data.txt                # Remove file
rm -r old_dir              # Remove directory and contents (-r = recursive)
```

### Wildcards

| Pattern | Matches |
|---------|---------|
| `*` | Any sequence of characters |
| `?` | Any single character |

```bash
cp *.fastq results/        # Copy all .fastq files
ls *.gz                    # List all gzipped files
rm temp_?.txt              # Remove temp_1.txt, temp_a.txt, etc.
```

### ⚠️ CAUTION with `rm -rf`

`rm -rf directory` **permanently deletes** a directory and everything in it. There is no undo. Double-check the path before pressing Enter. Never run `rm -rf /`—it can destroy your system.

---

## 4. Viewing File Contents

| Command | Use Case |
|---------|----------|
| `cat file` | Display entire file (small files) |
| `head file` | First 10 lines (default) |
| `head -n 15 file` | First 15 lines |
| `tail file` | Last 10 lines |
| `tail -n 20 file` | Last 20 lines |
| `less file` | Scroll through file interactively (q to quit) |
| `wc file` | Line, word, byte count |
| `wc -l file` | Line count only |
| `wc -w file` | Word count only |

**Bioinformatics tip:** Use `head` to peek at FASTA/FASTQ headers before processing.

---

## 5. Pipes and Redirects

### Streams

- **stdout** (1): Normal output
- **stderr** (2): Error messages
- **stdin** (0): Input

### Redirection

| Operator | Action | Example |
|----------|--------|---------|
| `>` | Overwrite file with output | `ls -l > file_list.txt` |
| `>>` | Append to file | `echo "done" >> log.txt` |
| `<` | Read input from file | `program < input.txt` |

### Pipe (`|`)

Sends stdout of one command as stdin to the next:

```bash
cat file.txt | wc -l           # Count lines
ls | grep fastq                # List only fastq files
sort data.txt | uniq | wc -l   # Sort, remove duplicates, count
```

**Bioinformatics example:** Count unique sequence IDs in a FASTA:

```bash
grep "^>" sequences.fasta | sort | uniq | wc -l
```

---

## 6. Pattern Matching with grep

`grep` searches for lines matching a pattern.

### Basic Usage

```bash
grep "pattern" file.txt       # Lines containing "pattern"
grep "error" server.log       # Find error lines
```

### Common Flags

| Flag | Effect |
|------|--------|
| `-i` | Case-insensitive |
| `-r` | Recursive (search in directories) |
| `-n` | Show line numbers |
| `-c` | Count matches per file |
| `-v` | Invert—lines that *don't* match |

```bash
grep -i error server.log                    # Case-insensitive
grep -rn "ATGC" sequences/                  # Recursive, with line numbers
grep -v "^#" config.txt                     # Exclude comment lines
```

### Simple Regex Patterns

- `^` = start of line
- `$` = end of line
- `.` = any character
- `*` = zero or more of preceding character

**Bioinformatics example:** Extract FASTA headers:

```bash
grep "^>" sequences.fasta
```

---

## 7. File Permissions

### Understanding `ls -l` Output

```bash
$ ls -l
-rwxr-xr-x 1 user group 1024 Jan 15 10:30 script.sh
drwxr-xr-x 2 user group 4096 Jan 15 10:30 data/
```

**First character:** File type

- `d` = directory  
- `-` = regular file  
- `l` = symbolic link  

**Next 9 characters:** Permissions (3 groups of 3)

- `rwx` = owner: read, write, execute
- `r-x` = group: read, execute
- `r-x` = others: read, execute

### chmod Basics

**Symbolic:**

```bash
chmod u+x script.sh     # Add execute for user (owner)
chmod g+w file.txt      # Add write for group
chmod o-r file.txt      # Remove read for others
```

**Numeric (common):**

- 4 = read, 2 = write, 1 = execute
- 7 = rwx, 6 = rw-, 5 = r-x

```bash
chmod 755 script.sh     # rwxr-xr-x (owner: all; others: read+execute)
chmod 644 file.txt      # rw-r--r--
```

**Making scripts executable:**

```bash
chmod +x myscript.sh
./myscript.sh
```

---

## 8. Environment Variables and PATH

### Displaying Variables

```bash
echo $PATH              # Show PATH (colon-separated directories)
echo $HOME              # Your home directory
```

### Setting Variables

```bash
export MY_VAR="value"   # Set for current session
```

### Configuration Files

- `~/.bashrc` — Runs for interactive non-login shells (typical terminal)
- `~/.bash_profile` — Runs for login shells

Add to PATH permanently:

```bash
echo 'export PATH="$PATH:$HOME/bin"' >> ~/.bashrc
source ~/.bashrc
```

### Why PATH Matters for Bioinformatics

Programs like `bwa`, `samtools`, or `bowtie2` must be in PATH to run by name. If you install tools in `~/bin`, add that to PATH.

---

## 9. Process Management

| Command | Purpose |
|---------|---------|
| `ps` | List processes (`ps aux` for detailed list) |
| `top` / `htop` | Interactive process viewer |
| `&` | Run command in background |
| `Ctrl+C` | Interrupt (kill) foreground process |
| `kill PID` | Terminate process by PID |
| `kill -9 PID` | Force kill (SIGKILL) |
| `nohup command &` | Run immune to hangups (survives logout) |

```bash
./long_analysis.sh &     # Run in background
kill 12345              # Terminate process 12345
nohup align_reads.sh > align.log 2>&1 &
```

---

## 10. Compression and Archives

### gzip / gunzip

```bash
gzip file.txt            # Creates file.txt.gz (removes original)
gunzip file.txt.gz       # Decompress
gzip -d file.txt.gz      # Same as gunzip
```

### tar (Archives)

| Operation | Flags | Example |
|-----------|-------|---------|
| Create | `-c` | `tar -czf archive.tar.gz dir/` |
| Extract | `-x` | `tar -xzf archive.tar.gz` |
| List | `-t` | `tar -tzf archive.tar.gz` |

**Common flags:** `-z` = gzip, `-f` = file

```bash
tar -czf project.tar.gz project/     # Create compressed archive
tar -xzf project.tar.gz               # Extract
```

### Common Bioinformatics Formats

- `.gz` — gzipped (FASTA, FASTQ, SAM, VCF)
- `.tar.gz` / `.tgz` — Tar + gzip (datasets, tools)

---

## 11. SSH and Remote Connections

### Connect to Remote Server

```bash
ssh scientist@lab.university.edu
```

### Copying Files (scp)

```bash
scp file.txt user@host:/remote/path/     # Upload
scp user@host:/remote/file.txt .         # Download
scp -r folder/ user@host:/remote/path/   # Recursive (directory)
```

### SSH Keys (Passwordless Login)

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
ssh-copy-id user@host    # Copy public key to server
```

---

## 12. Useful Extras

### awk One-Liners

```bash
awk '{print $1}' file.txt     # Print first column
awk -F',' '{print $2}' file   # Comma-separated, 2nd column
awk 'NR>1' file.txt           # Skip header line
```

### sed One-Liners

```bash
sed 's/old/new/' file.txt     # Substitute first occurrence per line
sed 's/old/new/g' file.txt   # Substitute all (global)
sed -i 's/old/new/g' file    # In-place edit (use with care)
```

### Downloading

```bash
curl -O https://example.com/file.fasta    # -O saves with original name
wget https://example.com/file.fasta       # Direct download
```

### Getting Help

```bash
man grep          # Manual page
grep --help       # Quick help
which samtools    # Where is this program?
```

---

## Practice Exercises

### Exercise 1: Navigation and Listing
Create a directory `bioinf_practice`, navigate into it, and list its contents with human-readable sizes. Then go back to the parent directory.

<details>
<summary>Solution</summary>

```bash
mkdir bioinf_practice
cd bioinf_practice
ls -lh
cd ..
```
</details>

### Exercise 2: File Operations and Redirects
Create a file `samples.txt` with three lines (sample1, sample2, sample3). Copy it to a subdirectory `backup`, then append a fourth line "sample4" to the original.

<details>
<summary>Solution</summary>

```bash
mkdir -p backup
echo -e "sample1\nsample2\nsample3" > samples.txt
cp samples.txt backup/
echo "sample4" >> samples.txt
```
</details>

### Exercise 3: Pipes and grep
Given a FASTA file `genes.fasta`, count how many sequences (headers) it contains. Hint: headers start with `>`.

<details>
<summary>Solution</summary>

```bash
grep -c "^>" genes.fasta
# OR
grep "^>" genes.fasta | wc -l
```
</details>

### Exercise 4: Viewing and Processing
Use `head` to view the first 5 lines of a large file, then use `wc -l` to count its total lines. Chain with a pipe if you can.

<details>
<summary>Solution</summary>

```bash
head -n 5 large_file.txt
wc -l large_file.txt
# For pipe (filtering): head doesn't change line count, so wc -l on full file
# A pipe example: cat large_file.txt | head -n 100 | wc -l  # counts 100
```
</details>

### Exercise 5: Permissions and Execution
Create a script `hello.sh` containing `#!/bin/bash` and `echo "Hello"`. Make it executable and run it.

<details>
<summary>Solution</summary>

```bash
echo -e '#!/bin/bash\necho "Hello"' > hello.sh
chmod +x hello.sh
./hello.sh
```
</details>

---

## Additional Resources

- [The Unix Shell (Software Carpentry)](https://swcarpentry.github.io/shell-novice/) — Free, interactive lessons
- [Linux Command Line Basics (Ubuntu)](https://ubuntu.com/tutorials/command-line-for-beginners) — Quick start
- [Explain Shell](https://explainshell.com/) — Paste a command to see what each part does
- [Bash Hackers Wiki](https://wiki.bash-hackers.org/) — Reference for advanced topics

---

*Good luck! Practice in a real terminal until these commands feel natural.*
