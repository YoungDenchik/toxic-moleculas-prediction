import type { LipinskiProperties as LipinskiPropertiesType } from "../types";

interface LipinskiPropertiesProps {
  properties: LipinskiPropertiesType;
  lipinskiPass: boolean;
}

const LIPINSKI_RULES = {
  molecular_weight: { name: "Molecular Weight", unit: "Da", threshold: "≤ 500" },
  logp: { name: "LogP", unit: "", threshold: "≤ 5" },
  hbd: { name: "H-Bond Donors", unit: "", threshold: "≤ 5" },
  hba: { name: "H-Bond Acceptors", unit: "", threshold: "≤ 10" },
};

export function LipinskiProperties({ properties, lipinskiPass }: LipinskiPropertiesProps) {
  const propertyKeys = Object.keys(LIPINSKI_RULES) as (keyof typeof LIPINSKI_RULES)[];

  return (
    <div style={{ border: "1px solid #e2e8f0", borderRadius: 16, padding: 16 }}>
      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", marginBottom: 12 }}>
        <div style={{ fontWeight: 800, fontSize: 15 }}>Lipinski Rule of 5</div>
        <span
          style={{
            padding: "4px 10px",
            borderRadius: 999,
            fontSize: 12,
            fontWeight: 700,
            background: lipinskiPass ? "#dcfce7" : "#fee2e2",
            color: lipinskiPass ? "#166534" : "#991b1b",
          }}
        >
          {lipinskiPass ? "PASS" : "FAIL"}
        </span>
      </div>

      <div style={{ display: "grid", gap: 8 }}>
        {propertyKeys.map((key) => {
          const prop = properties[key];
          const rule = LIPINSKI_RULES[key];
          const passColor = prop.pass_rule ? "#16a34a" : "#dc2626";
          const passBg = prop.pass_rule ? "#f0fdf4" : "#fef2f2";

          return (
            <div
              key={key}
              style={{
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
                padding: "8px 12px",
                borderRadius: 8,
                background: passBg,
                border: `1px solid ${prop.pass_rule ? "#bbf7d0" : "#fecaca"}`,
              }}
            >
              <div style={{ display: "flex", flexDirection: "column" }}>
                <span style={{ fontWeight: 600, fontSize: 13, color: "#0f172a" }}>{rule.name}</span>
                <span style={{ fontSize: 11, color: "#64748b" }}>Rule: {rule.threshold}</span>
              </div>
              <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                <span
                  style={{
                    fontWeight: 700,
                    fontSize: 14,
                    color: passColor,
                  }}
                >
                  {typeof prop.value === "number"
                    ? prop.value % 1 === 0
                      ? prop.value
                      : prop.value.toFixed(2)
                    : prop.value}
                  {rule.unit && <span style={{ fontSize: 11, fontWeight: 400 }}> {rule.unit}</span>}
                </span>
                <span style={{ fontSize: 16 }}>{prop.pass_rule ? "✓" : "✗"}</span>
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}
