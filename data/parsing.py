import json
from collections import defaultdict
    
def offline(instance):
    data = {}
    with open(instance+"/offline.json") as f:
        offline = json.load(f)
    data["o"] = "O"
    while data["o"] in offline["shifts_length"].keys():
        data["o"] += "ff"
    data["Sigma"] = set(offline["shifts_length"].keys()).union({data["o"]})
    data["tau"] = {**offline["shifts_length"],**{data["o"]:0}}
    data["G"] = set(offline["contracts"].keys())
    data["F"] = {g:{s:set(contract["rotation"].get(s,[])) for s in offline["shifts_length"].keys()} for g,contract in offline["contracts"].items()}
    data["gamma_max"] = {g:contract["stretchs"]["max_on"] for g,contract in offline["contracts"].items()}
    data["gamma_min"] = {g:contract["stretchs"]["min_on"] for g,contract in offline["contracts"].items()}
    data["delta_min"] = {g:contract["stretchs"]["min_off"] for g,contract in offline["contracts"].items()}
    data["nu_max"] = {g:contract["weekends"] for g,contract in offline["contracts"].items()}
    data["tau_max"] = {g:contract["max_workload"] for g,contract in offline["contracts"].items()}    
    data["rho_max"] = {g:{s:contract["max_shifts"].get(s,float("inf")) for s in offline["shifts_length"].keys()} for g,contract in offline["contracts"].items()}
    return data

def online(instance,version = "personalized"):
    data = {}
    with open(instance+"/online.json") as f:
        online = json.load(f)
    data["E"] = defaultdict(set)
    data["T"] = online["horizon"]
    data["W"] = data["T"]//7
    data["d"] = {d["day"]:d["required"] for d in online["demand"]}
    data["u"] = {d["day"]:{s:online["costs"]["undercover"] for s in d["required"]} for d in online["demand"]}
    data["v"] = {d["day"]:{s:online["costs"]["overcover"] for s in d["required"]} for d in online["demand"]}
    data["O"] = {}
    data["p"] = {}
    data["q"] = {}
    for k,v in online["employees"].items():
        data["E"][v["contract"]].add(k)
        data["p"][k] = {t:{s:0 for s in d.keys()} for t,d in data["d"].items()}
        data["q"][k] = {t:{s:0 for s in d.keys()} for t,d in data["d"].items()}
        if version == "personalized":
            data["O"][k] = set(v["days_off"])
            for r in v["request_on"]:
                data["p"][k][r["day"]][r["shift"]] = r["weight"]
            for r in v["request_off"]:
                data["q"][k][r["day"]][r["shift"]] = r["weight"]
        elif version == "anonymous":
            data["O"][k] = set()
        else:
            raise Exception("Version has to be either 'personalized' or 'anonymous'") 
    data["E"] = dict(data["E"])
    return data