# Course 7: Glossary

Short definitions for terms used in Part 1 (core inference) and Part 2 (biology application).

| Term | Definition |
|------|------------|
| **Batching** | Grouping multiple inputs (e.g. requests, sequences) to be processed together on hardware, improving throughput and GPU utilization. |
| **Inference** | Running a trained model to produce outputs (predictions, embeddings, etc.); contrasted with **training**. |
| **Training** | The phase where model parameters are learned from data; usually more compute- and memory-intensive than inference. |
| **GPU saturation** | Using the GPU’s compute and memory bandwidth heavily so that it is not idle; a goal of efficient inference. |
| **MIG (Multi-Instance GPU)** | NVIDIA feature to partition a single GPU into smaller instances so multiple jobs or users can share the GPU with isolation. |
| **Sharding** | Splitting a model or workload across multiple devices or nodes (e.g. different GPUs) to fit memory or increase throughput. |
| **Parallelization** | Using multiple compute units (cores, GPUs, nodes) at once to speed up or scale inference. |
| **Memory management (inference)** | Controlling how activations, weights, and KV caches are stored and moved so that inference fits in available memory and runs efficiently. |
| **Concurrent user sessions** | Multiple users or requests being served at the same time; requires batching, scheduling, and often multi-GPU or MIG. |
| **KV cache** | Stored key-value pairs from previous tokens in autoregressive decoding, to avoid recomputation; a major consumer of memory in LLM inference. |
| **Throughput** | Amount of work (e.g. tokens or requests) completed per unit time; often the main metric for inference systems. |
| **Latency** | Time to complete a single request or token; important for interactive use. |
| **Gene-language model** | A model trained on gene or transcript sequences (e.g. gene expression, sequences); inference concepts apply when serving such models. |
| **Cell embedding model** | A model that maps cells (e.g. from single-cell data) to vectors; batching and GPU use matter when embedding large atlases. |
