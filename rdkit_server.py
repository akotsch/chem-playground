from http.server import BaseHTTPRequestHandler, HTTPServer
import json
import base64
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Draw


# ===============================
# CHEMISTRY LOGIC (OUTSIDE CLASS)
# ===============================

def suggest_arrows(reactant_smiles, reagent_smiles):
    arrows = []

    mol = Chem.MolFromSmiles(reactant_smiles)
    reagent = Chem.MolFromSmiles(reagent_smiles)

    if not mol or not reagent:
        return arrows

    # Rule: Alkene + Br2
    alkene = Chem.MolFromSmarts("C=C")
    bromine = Chem.MolFromSmarts("BrBr")

    if mol.HasSubstructMatch(alkene) and reagent.HasSubstructMatch(bromine):
        a1, a2 = mol.GetSubstructMatch(alkene)

        arrows.append({
            "start_atom": a1,
            "end_atom": a2,
            "type": "pi_attack"
        })

    return arrows


# ===============================
# HTTP SERVER
# ===============================

class Handler(BaseHTTPRequestHandler):

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")

    def do_OPTIONS(self):
        self.send_response(200)
        self._cors()
        self.end_headers()

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        data = json.loads(self.rfile.read(length))

        # ---------- RENDER ----------
        if self.path == "/render":
            smiles = data.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                self.send_response(400)
                self._cors()
                self.end_headers()
                self.wfile.write(b"Invalid SMILES")
                return

            img = Draw.MolToImage(mol, size=(350, 250))
            buf = BytesIO()
            img.save(buf, format="PNG")

            payload = json.dumps({
                "image": base64.b64encode(buf.getvalue()).decode()
            }).encode()

            self.send_response(200)
            self._cors()
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return

        # ---------- SUGGEST (LEVEL 1) ----------
        if self.path == "/suggest":
            reactant = data.get("reactant", "")
            reagent = data.get("reagent", "")

            arrows = suggest_arrows(reactant, reagent)

            payload = json.dumps({
                "arrows": arrows
            }).encode()

            self.send_response(200)
            self._cors()
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return

        # ---------- FALLBACK ----------
        self.send_response(404)
        self.end_headers()


print("RDKit server running at http://localhost:8000")
HTTPServer(("localhost", 8000), Handler).serve_forever()
