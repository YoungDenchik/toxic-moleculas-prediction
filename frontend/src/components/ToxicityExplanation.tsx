import type { Explainability, Prediction } from "../types";

interface ToxicityExplanationProps {
  prediction: Prediction;
  explainability: Explainability | null;
}

const TOXICITY_ICONS: Record<string, string> = {
  aromatic: "🔬",
  amine: "⚗️",
  quinone: "⚠️",
  lipophilicity: "💧",
  molecular: "⚖️",
  hydrogen: "🔗",
  default: "📊",
};

function getIconForFactor(factor: string): string {
  const lowerFactor = factor.toLowerCase();
  if (lowerFactor.includes("aromatic")) return TOXICITY_ICONS.aromatic;
  if (lowerFactor.includes("amine")) return TOXICITY_ICONS.amine;
  if (lowerFactor.includes("quinone")) return TOXICITY_ICONS.quinone;
  if (lowerFactor.includes("lipophilicity") || lowerFactor.includes("logp")) return TOXICITY_ICONS.lipophilicity;
  if (lowerFactor.includes("molecular weight") || lowerFactor.includes("clearance")) return TOXICITY_ICONS.molecular;
  if (lowerFactor.includes("hydrogen") || lowerFactor.includes("h-bond")) return TOXICITY_ICONS.hydrogen;
  return TOXICITY_ICONS.default;
}

export function ToxicityExplanation({ prediction, explainability }: ToxicityExplanationProps) {
  const isToxic = prediction.label === "toxic";
  const probability = prediction.probability;

  const riskLevel =
    probability >= 0.8 ? "Very High" : probability >= 0.6 ? "High" : probability >= 0.4 ? "Moderate" : "Low";

  const riskColor =
    probability >= 0.8
      ? "#991b1b"
      : probability >= 0.6
        ? "#dc2626"
        : probability >= 0.4
          ? "#ea580c"
          : "#16a34a";

  const riskBg =
    probability >= 0.8
      ? "#fef2f2"
      : probability >= 0.6
        ? "#fff1f2"
        : probability >= 0.4
          ? "#fff7ed"
          : "#f0fdf4";

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
      <div style={{ fontWeight: 800, fontSize: 15, marginBottom: 12 }}>Cardiotoxicity Prediction</div>

      {/* Main prediction result */}
      <div
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          padding: "12px 16px",
          borderRadius: 12,
          background: isToxic ? "#fef2f2" : "#f0fdf4",
          border: `2px solid ${isToxic ? "#fecaca" : "#bbf7d0"}`,
          marginBottom: 16,
        }}
      >
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <span style={{ fontSize: 32 }}>{isToxic ? "⚠️" : "✅"}</span>
          <div>
            <div
              style={{
                fontWeight: 800,
                fontSize: 18,
                color: isToxic ? "#991b1b" : "#166534",
              }}
            >
              {isToxic ? "TOXIC" : "NON-TOXIC"}
            </div>
            <div style={{ fontSize: 12, color: "#64748b" }}>hERG channel interaction risk</div>
          </div>
        </div>
        <div style={{ textAlign: "right" }}>
          <div style={{ fontWeight: 700, fontSize: 24, color: riskColor }}>{(probability * 100).toFixed(1)}%</div>
          <div
            style={{
              fontSize: 11,
              fontWeight: 600,
              color: riskColor,
              background: riskBg,
              padding: "2px 8px",
              borderRadius: 999,
            }}
          >
            {riskLevel} Risk
          </div>
        </div>
      </div>

      {/* Probability bar */}
      <div style={{ marginBottom: 16 }}>
        <div style={{ fontSize: 12, color: "#64748b", marginBottom: 4 }}>Toxicity Probability</div>
        <div
          style={{
            height: 8,
            background: "#e2e8f0",
            borderRadius: 4,
            overflow: "hidden",
          }}
        >
          <div
            style={{
              height: "100%",
              width: `${probability * 100}%`,
              background: `linear-gradient(90deg, #22c55e 0%, #eab308 50%, #ef4444 100%)`,
              borderRadius: 4,
              transition: "width 0.3s ease",
            }}
          />
        </div>
        <div style={{ display: "flex", justifyContent: "space-between", fontSize: 10, color: "#94a3b8", marginTop: 2 }}>
          <span>Safe</span>
          <span>Toxic</span>
        </div>
      </div>

      {/* Explanation factors */}
      {isToxic && explainability && explainability.factors.length > 0 && (
        <div>
          <div
            style={{
              fontWeight: 700,
              fontSize: 13,
              marginBottom: 8,
              color: "#0f172a",
              display: "flex",
              alignItems: "center",
              gap: 6,
            }}
          >
            <span>🔍</span> Why might this molecule be toxic?
          </div>
          <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
            {explainability.factors.map((factor, idx) => (
              <div
                key={idx}
                style={{
                  display: "flex",
                  alignItems: "flex-start",
                  gap: 10,
                  padding: "10px 12px",
                  background: "#fffbeb",
                  border: "1px solid #fde68a",
                  borderRadius: 8,
                }}
              >
                <span style={{ fontSize: 18, flexShrink: 0 }}>{getIconForFactor(factor)}</span>
                <span style={{ fontSize: 13, color: "#78350f", lineHeight: 1.4 }}>{factor}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Safe molecule message */}
      {!isToxic && (
        <div
          style={{
            padding: "12px 16px",
            background: "#f0fdf4",
            border: "1px solid #bbf7d0",
            borderRadius: 8,
          }}
        >
          <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
            <span style={{ fontSize: 18 }}>🛡️</span>
            <span style={{ fontSize: 13, color: "#166534", fontWeight: 500 }}>
              This molecule shows low hERG channel binding risk based on its structural features and physicochemical
              properties.
            </span>
          </div>
        </div>
      )}
    </div>
  );
}
