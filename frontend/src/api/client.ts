import type { PredictRequest, PredictResponse, UmapPoint } from "../types";

const BASE = (import.meta.env.VITE_API_BASE_URL as string) || "";

async function fetchJSON<T>(url: string, options?: RequestInit): Promise<T> {
  const res = await fetch(url, {
    ...options,
    headers: { "Content-Type": "application/json", ...(options?.headers || {}) },
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
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
