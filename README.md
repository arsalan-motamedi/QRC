# Quantumness in Reservoir Computing

Welcome to the GitHub repository dedicated to the paper entitled "Correlations Between Quantumness and Learning Performance in Reservoir Computing with a Single Oscillator." This paper was a collaborative effort between myself, [Dr. Hadi Zadeh-Haghighi](https://contacts.ucalgary.ca/info/phas/profiles/1-9226636), and [Prof. Christoph Simon](https://science.ucalgary.ca/physics-astronomy/contacts/christoph-simon) from the University of Calgary.

We have made the decision to share the code used in our research with you. Please feel free to reach out to us if you have any questions or comments regarding the contents of this repository.

To give you a brief overview of our work, we investigate the impact of quantumness in a reservoir computing task. Our reservoir consists of a simple Kerr oscillator, and we focus on predicting Mackey-Glass time-series. Our results suggest that quantumness could be a valuable resource in the learning process. For more details, please refer to [the paper](https://arxiv.org/abs/2304.03462).

## Learning Performance
In Figure 2 folder, you can find the main learnning codes, for different tasks, such as learning a Mackey-Glass or a Rossler dynamics. Below is the training result for the Mackey-Glass.

- Mackey-Glass dynamics:

<img src="https://user-images.githubusercontent.com/94669474/229359616-d795df4d-7195-4aa5-a757-38171729136b.jpg" width=25% height=25%> <img src="https://user-images.githubusercontent.com/94669474/229359717-9fdce73e-48f7-4f1f-a684-5b1ce0fd018a.jpg" width=25% height=25%>

- Rossler dynamics:

<img src="https://user-images.githubusercontent.com/94669474/229359833-83fc06a8-a44d-42fe-9f57-646a7a02e8f9.jpg" width=25% height=25%> <img src="https://user-images.githubusercontent.com/94669474/229359846-b97054de-30df-4f7c-82a0-14cc98102be8.jpg" width=25% height=25%>

- Noisy periodic functions:

<img src="https://user-images.githubusercontent.com/94669474/229360384-d3267cb2-dc40-40e9-8990-a1b5d8d5f643.jpg" width=40% height=40%>

## The effect of quantumness

We use a set of 35 random states and 30 different hyperparameters in the training of the reservoir, and examine the interplay between the quantumness of the reservoir and the learning performance. Using figures below, we deduce that quantumness is indeed a game-changer! We first note that a more quantum state 'can' improve the performance. We furthermore, identify conditions on our hyperparameters that give rise to the optimal performance, and further, we discuss that these results match intuition. Please read [our paper](https://arxiv.org/abs/2304.03462) for a detailed discussion.

<img src="https://github.com/arsalan-motamedi/QRC/assets/94669474/0cb6421f-d965-4e4d-ab0e-f674dc9b29ca" width=50% height=50%>


