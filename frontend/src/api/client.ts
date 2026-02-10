import type {
  PredictRequest,
  PredictResponse,
  UmapPoint,
  AnalyzeResponse,
  UmapResponse,
  LipinskiResponse,
  ExplainResponse,
} from "../types";

const BASE = (import.meta.env.VITE_API_BASE_URL as string) || "";

// Debug: check what URL is being used
console.log("[API] VITE_API_BASE_URL:", import.meta.env.VITE_API_BASE_URL);
console.log("[API] BASE URL:", BASE || "(empty - using same origin)");

async function fetchJSON<T>(url: string, options?: RequestInit): Promise<T> {
  const res = await fetch(url, {
    ...options,
    headers: { "Content-Type": "application/json", ...(options?.headers || {}) },
  });
  if (!res.ok) {
    const errorData = await res.json().catch(() => ({}));
    throw new Error(errorData.detail || `HTTP ${res.status}`);
  }
  return (await res.json()) as T;
}

export async function apiGetUmap(): Promise<UmapPoint[]> {
  return fetchJSON<UmapPoint[]>(`${BASE}/umap`);
}

export async function apiPredict(body: PredictRequest): Promise<PredictResponse> {
  return fetchJSON<PredictResponse>(`${BASE}/predict`, {
    method: "POST",
    body: JSON.stringify(body),
  });
}

export async function apiAnalyze(smiles: string): Promise<AnalyzeResponse> {
  return fetchJSON<AnalyzeResponse>(`${BASE}/analyze`, {
    method: "POST",
    body: JSON.stringify({ smiles }),
  });
}

export async function apiGetUmapProjection(smiles: string): Promise<UmapResponse> {
  return fetchJSON<UmapResponse>(`${BASE}/chemical-space/umap`, {
    method: "POST",
    body: JSON.stringify({ smiles }),
  });
}

export async function apiGetLipinski(smiles: string): Promise<LipinskiResponse> {
  return fetchJSON<LipinskiResponse>(`${BASE}/properties/lipinski`, {
    method: "POST",
    body: JSON.stringify({ smiles }),
  });
}

export async function apiExplainToxicity(smiles: string): Promise<ExplainResponse> {
  return fetchJSON<ExplainResponse>(`${BASE}/explain/toxicity`, {
    method: "POST",
    body: JSON.stringify({ smiles }),
  });
}
