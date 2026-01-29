import { useState } from "react";

export function PredictPanel() {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<null | { risk: number; label: "LOW" | "HIGH" }>(null);
  const [error, setError] = useState<string | null>(null);

  async function onPredict() {
    setError(null);
    setResult(null);

    const s = smiles.trim();
    if (!s) {
      setError("Введи SMILES (рядок молекули).");
      return;
    }

    // Поки бек не підключений — робимо фейковий результат
    setLoading(true);
    try {
      await new Promise((r) => setTimeout(r, 600));
      const fake = Math.random(); // 0..1
      setResult({
        risk: fake,
        label: fake >= 0.5 ? "HIGH" : "LOW",
      });
    } catch {
      setError("Не вдалося зробити прогноз.");
    } finally {
      setLoading(false);
    }
  }

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 14, overflow: "hidden" }}>
      <div style={{ fontWeight: 800, marginBottom: 10 }}>Predict</div>

      <label style={{ display: "block", fontSize: 12, color: "#64748b", marginBottom: 6 }}>
        SMILES
      </label>

      <textarea
        value={smiles}
        onChange={(e) => setSmiles(e.target.value)}
        placeholder='Напр: "CCO" або встав свій SMILES...'
        rows={4}
        style={{
          width: "100%",
          boxSizing: "border-box",
          resize: "vertical",
          border: "1px solid #cbd5e1",
          borderRadius: 12,
          padding: 10,
          fontFamily: "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace",
          fontSize: 13,
          outline: "none",
        }}
      />

      <div style={{ display: "flex", gap: 10, marginTop: 10, alignItems: "center" }}>
        <button
          onClick={onPredict}
          disabled={loading}
          style={{
            background: "#1d4ed8",
            color: "white",
            padding: "10px 14px",
            borderRadius: 12,
            border: "none",
            fontWeight: 700,
            cursor: loading ? "not-allowed" : "pointer",
          }}
        >
          {loading ? "Predicting..." : "Predict"}
        </button>

        <button
          onClick={() => {
            setSmiles("");
            setResult(null);
            setError(null);
          }}
          style={{
            background: "white",
            color: "#0f172a",
            padding: "10px 14px",
            borderRadius: 12,
            border: "1px solid #cbd5e1",
            fontWeight: 700,
            cursor: "pointer",
          }}
        >
          Clear
        </button>
      </div>

      {error && (
        <div style={{ marginTop: 10, color: "#b91c1c", fontSize: 13 }}>
          {error}
        </div>
      )}

      {result && (
        <div style={{ marginTop: 12, borderTop: "1px solid #e2e8f0", paddingTop: 12 }}>
          <div style={{ fontWeight: 800, marginBottom: 6 }}>Result</div>
          <div style={{ display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap" }}>
            <span
              style={{
                padding: "6px 10px",
                borderRadius: 999,
                border: "1px solid #cbd5e1",
                background: "#f8fafc",
                fontSize: 12,
                fontWeight: 800,
              }}
            >
              {result.label}
            </span>

            <span style={{ color: "#334155", fontSize: 13 }}>
              risk score: <b>{result.risk.toFixed(3)}</b>
            </span>
          </div>

          <div style={{ marginTop: 8, fontSize: 12, color: "#64748b" }}>
            (поки це mock. Коли бек буде готовий — замінимо на справжній API.)
          </div>
        </div>
      )}
    </div>
  );
}
