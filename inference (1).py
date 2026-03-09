import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem

# ------------------------
# Load trained models
# ------------------------
rf_tc = joblib.load("rf_tc.pkl")
rf_tg = joblib.load("rf_tg.pkl")
rf_density = joblib.load("rf_density.pkl")
rf_ffv = joblib.load("rf_ffv.pkl")
rf_rg = joblib.load("rf_rg.pkl")

N_BITS = 2048

# ------------------------
# Fingerprint function
# ------------------------
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

N_BITS = 2048

def smiles_to_fingerprint(smiles, n_bits=N_BITS):

    try:
        # Fix polymer wildcard atoms
        smiles = smiles.replace("*", "C")

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print("RDKit failed:", smiles)
            return None

        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=2,
            nBits=n_bits
        )

        return np.array(fp).reshape(1, -1)

    except Exception as e:
        print("Fingerprint error:", e)
        return None

# ------------------------
# Single polymer prediction
# ------------------------
def predict_properties(smiles):

    fp = smiles_to_fingerprint(smiles)

    if fp is None:
        return None   # IMPORTANT change

    return {
        "Tc": float(rf_tc.predict(fp)[0]),
        "Tg": float(rf_tg.predict(fp)[0]),
        "Density": float(rf_density.predict(fp)[0]),
        "FFV": float(rf_ffv.predict(fp)[0]),
        "Rg": float(rf_rg.predict(fp)[0])
    }

# ------------------------
# Blend fingerprint
# ------------------------
def blend_fingerprint(smiles_a, smiles_b, w=0.5):
    fp_a = smiles_to_fingerprint(smiles_a)
    fp_b = smiles_to_fingerprint(smiles_b)

    if fp_a is None or fp_b is None:
        return None

    return w * fp_a + (1 - w) * fp_b

# ------------------------
# Blend prediction
# ------------------------
def predict_blend_properties(smiles_a, smiles_b, w=0.5):
    fp = blend_fingerprint(smiles_a, smiles_b, w)
    if fp is None:
        return {"error": "Invalid SMILES"}

    return {
        "Tc": float(rf_tc.predict(fp)[0]),
        "Tg": float(rf_tg.predict(fp)[0]),
        "Density": float(rf_density.predict(fp)[0]),
        "FFV": float(rf_ffv.predict(fp)[0]),
        "Rg": float(rf_rg.predict(fp)[0]),
    }
def check_sustainability(pred):
    """
    Rule-based sustainability check using predicted properties
    """

    score = 0

    # Rules (you can adjust thresholds)
    if pred["Tg"] >= 60:        # thermal stability
        score += 1

    if pred["Density"] <= 1.2: # lightweight material
        score += 1

    if pred["FFV"] >= 0.35:    # good free volume
        score += 1

    if pred["Rg"] <= 20:       # compact structure
        score += 1

    if score >= 3:
        return "Sustainable"
    else:
        return "Not Sustainable"

def predict_blend_with_sustainability(smiles_a, smiles_b, w=0.5):
    pred = predict_blend_properties(smiles_a, smiles_b, w)

    # 🔴 FIX: handle invalid predictions
    if "error" in pred:
        pred["Sustainability"] = "Invalid"
        return pred

    status = check_sustainability(pred)
    pred["Sustainability"] = status
    return pred




 


