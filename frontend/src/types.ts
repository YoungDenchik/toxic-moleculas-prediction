export type UmapPoint = {
  id: string;
  x: number;
  y: number;
  label?: number;
  smiles?: string;
  risk?: "Low" | "Medium" | "High";
};

export type PredictRequest = { smiles: string };

export type PropertyValue = {
  value: number;
  pass_rule: boolean;
};

export type LipinskiProperties = {
  molecular_weight: PropertyValue;
  logp: PropertyValue;
  hbd: PropertyValue;
  hba: PropertyValue;
};

export type Prediction = {
  label: "toxic" | "non-toxic";
  probability: number;
};

export type ChemicalSpace = {
  x: number;
  y: number;
};

export type Explainability = {
  probability: number;
  factors: string[];
};

export type AnalyzeResponse = {
  canonical_smiles: string;
  prediction: Prediction;
  chemical_space: ChemicalSpace;
  properties: LipinskiProperties;
  lipinski_pass: boolean;
  explainability: Explainability | null;
};

export type PredictResponse = {
  canonical_smiles: string;
  prediction: Prediction;
};

export type UmapResponse = {
  method: string;
  canonical_smiles: string;
  prediction: Prediction;
  point: ChemicalSpace;
};

export type LipinskiResponse = {
  canonical_smiles: string;
  molecular_weight: PropertyValue;
  logp: PropertyValue;
  hbd: PropertyValue;
  hba: PropertyValue;
  lipinski_pass: boolean;
};

export type ExplainResponse = {
  probability: number;
  factors: string[];
};
