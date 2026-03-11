# Social Network Analysis : ISIS Network

This project aims at applying SNA technique to explore the relational, structural and emotional dynamics of a network that consists in supporters of ISIS on Twitter.\\
The analysis is based on a dataset composed by over 17000 tweets (Kaggle: *"How ISIS Uses Twitter"*), published between the Paris terroristic attack in 2015 and May 2016.

## 🎯 Goals and Key Results
The analysis reveals a sparse network that is strongly structured around few central users, characterized by a strong local influence and emotive component.

* **Key Nodes & Centrality:** * The most influent and mentioned account have results to be `RamiAlLolah`, `Nidalgazaui` e `WarReporter1`.
  * The analysis of *Betweenness Centrality* evidenced that figures as `Nidalgazaui` and `Uncle_SamCoco` happens to be crucila for the flow of information through different parts of the network.
* **Temporal Analysis:** It has been observed a strong tendency to **triadic closure**. Users have more probabiliy to create new links if they already have a common link, replying the existence of standard social dynamics even in radical context.
* **Community Detection:** thanks to the **Louvain** algorithm, there have been detected 8 principal communities. Between those, 3 "core-communities" dominates the flow of information, showing an high internal density and strong assortativity (users communicate prevalently inside their community).
* **Text & Emotion Analysis:** thanks to the use of '*NRC Emotion Lexicon*, it appears that messages are dominated by feeling of **Fear, Anger and Trust**. The most used hashtags (#isis, #syria, #aleppo) and the key words reflect a narrative focused on war events and identity legitimacy.

## 📊 Explore the Project

👉 **[Visualize the interactive Notebook on Kaggle](https://www.kaggle.com/code/pizzox/isis-social-network-analysis)** *(Suggested to explore code and graphics)*

📄 **[Read the Report in PDF format (in Italian)](docs/ISIS_SNA_Report.pdf)** 

---

## 📁 Structure of the Repository
* `docs/`: Contains the final report (PDF) in italian with the complete analysis.
* `src/`: Contains the script for the elaboration of the network the computation of centrality and the NLP.
* `data/`: Contains the dataset used in the analysis, the original dataset is provided by *Fifth Tribe* and can be downloaded from [Kaggle: How ISIS uses Twitter](https://www.kaggle.com/datasets/fifthtribe/how-isis-uses-twitter).