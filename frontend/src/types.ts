export type UmapPoint = {
  id: string;
  x: number;
  y: number;
  label?: number;
  smiles?: string;
  risk?: "Low" | "Medium" | "High";
};

export type PredictRequest = { smiles: string };

export type PredictResponse = {
  risk: "Low" | "Medium" | "High";
  probability: number;
  explanation: string[];
  umap?: { x: number; y: number };
};
