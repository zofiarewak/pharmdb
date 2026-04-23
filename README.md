# PharmDB – Drug Database & Analysis System
### Overview
PharmDB is a Python-based system for managing and analyzing a database of drugs.
It was developed as a final project for the Algorithms and Data Structures course (University of Warsaw, 2025).
### Features
Add drugs with:
- therapeutic indications (with efficacy score)
- side effects (severity, frequency)
- substitute relationships
### Fast analytics
#### O(1) queries for:
- number of indications above a threshold
- risk score of a drug
- worst side effect
- number of alternative drugs
### Run with
```
python pharmdb.py
```
### Example usage: 
```
db = PharmDB()

id1 = db.add_drug("Ibuprofen", [("pain", 8)], [], [("nausea", 2, 3)])
id2 = db.add_drug("Paracetamol", [("pain", 7)], [id1], [])

print(db.find_best_drug_for_indication("pain"))
```
