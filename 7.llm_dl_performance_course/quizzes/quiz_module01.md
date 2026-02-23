# Quiz: Part 1 — Core Inference (Modules 1–5)

**Target: ≥70% to pass. Mix of MC and short answer.**

---

## Multiple choice

### Q1
Why is efficient inference often cited as a major barrier to LLMs in production?

- A) Training is more expensive than inference  
- B) Throughput, memory, and latency constraints make serving at scale hard  
- C) GPUs are only used for training  
- D) Batching is not possible for LLMs  

**Your answer:** _____

---

### Q2
What does batching improve in inference?

- A) Only model accuracy  
- B) Throughput and GPU utilization by processing multiple inputs together  
- C) Only latency  
- D) Only memory use  

**Your answer:** _____

---

### Q3
What is MIG (Multi-Instance GPU)?

- A) A type of CPU  
- B) Partitioning one GPU into smaller instances for multiple jobs with isolation  
- C) A memory type  
- D) A batching strategy  

**Your answer:** _____

---

### Q4
GPU saturation refers to:

- A) Overheating the GPU  
- B) Keeping the GPU busy (high utilization of compute and memory bandwidth)  
- C) Using only one GPU  
- D) Turning off the GPU  

**Your answer:** _____

---

### Q5
Sharding in inference usually means:

- A) Caching results only  
- B) Splitting model or workload across multiple devices (e.g. GPUs)  
- C) Batching requests only  
- D) Using a single GPU  

**Your answer:** _____

---

## Short answer

### Q6
In one or two sentences, what is the main difference between inference and training in terms of compute and memory?

**Your answer:** _____

---

### Q7
Give one reason why managing memory is critical for LLM inference (e.g. what consumes memory?).

**Your answer:** _____

---

### Q8
When might you use MIG instead of giving each user a full GPU?

**Your answer:** _____

---

**Score: _____ / 8. Pass if ≥70% (e.g. 6/8 or better).**
