import requests
import pandas as pd
from core import config as cfg
import time

# Ova skripta vrši preuzimanje podataka iz BioGRID baze podataka, na osnovu kojih pravimo "PPI" mrežu
request_url = cfg.BASE_URL + "/interactions"
# spisak gena za koje dohvatamo podatke
geneList = [
    "AIP", "LEP", "CHGA", "CDKN2A", "CAV1", "CTNNB1", "ABCB1", "CASR", "ALB", "G6PD",
    "FLNA", "ESR1", "KL", "ACTB", "MEN1", "TP53", "KDR", "TICAM2", "EGFR", "CDKN2B",
    "EDNRB", "VDR", "LGALS3", "LRP5", "SYT1", "ENPEP", "MAPK1", "CCND1", "TNFSF10",
    "PTH", "PRKAR1A", "SLC12A3", "POTEF"
]

# eksperimente ovog tipa izbacujemo
evidenceList = ["Co-localization", "Genetic interference", "Synthetic Rescue", "Synthetic Growth Defect", "Synthetic Lethality"]

params = {
    "accesskey": cfg.ACCESS_KEY,
    "format": "tab3",   
    "geneList": "|".join(geneList),
    "interSpeciesExcluded": "true",
    "searchNames": "true",
    "evidenceList": "|".join(evidenceList), 
    "includeInteractors": "true",
    "includeInteractorInteractions": "true",
    "taxId": 9606,
    "includeEvidence": "false",
    "includeHeader": "true",
    "selfInteractionsExcluded": "true",
    "max": 10000
}

start = 0
all_interactions = ""
interaction_count = 0
batch_size = 100000
file_count = 0
retry_limit = 3  # koliko puta da se pokuša ukoliko dodje do greške
retry_delay = 10  # koliko sekundi da se čeka ukoliko dodje do greške da se ponovo pokuša sa preuzimanjem

while True:
    print(f"Fetching data starting at {start}")
    params["start"] = start
    success = False
    retries = 0

    while not success and retries < retry_limit:
        try:
            r = requests.get(request_url, params=params)
            r.raise_for_status()
            interactions = r.text
            params['includeHeader'] = "false"

            if not interactions.strip():
                break

            all_interactions += interactions
            interaction_count += interactions.count('\n')  

            if interaction_count >= batch_size-1:
                output_file = f"complete_network_part_{file_count}.txt"
                with open(output_file, "w") as file:
                    file.write(all_interactions)
                all_interactions = ""
                interaction_count = 0
                file_count += 1
                params['includeHeader'] = "true"

            start += params["max"]
            params['includeHeader'] = "false"
            success = True
        except requests.exceptions.HTTPError as e:
            print(f"HTTP error occurred: {e}. Retrying in {retry_delay} seconds...")
            retries += 1
            time.sleep(retry_delay)
        except Exception as e:
            print(f"An error occurred: {e}. Retrying in {retry_delay} seconds...")
            retries += 1
            time.sleep(retry_delay)
    
    if not success:
        print("Max retries reached. Exiting...")
        break

if all_interactions:
    output_file = f"complete_network_part_{file_count}.txt"
    with open(output_file, "w") as file:
        file.write(all_interactions)

print(f"Complete network saved to files with prefix 'complete_network_part_'")
