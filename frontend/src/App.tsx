import { useEffect, useState } from "react";
import type { ReactNode } from "react";
import type { CSSProperties } from "react";

import { apiGetUmap } from "./api/client";
import { MoleculeAnalyzer } from "./components/MoleculeAnalyzer";
import type { UmapPoint } from "./types";
import { mockUmap } from "./mocks/umap";

export default function App() {
  const [points, setPoints] = useState<UmapPoint[]>(mockUmap);
  const [loading, setLoading] = useState(true);
  const [usingMockData, setUsingMockData] = useState(false);

  useEffect(() => {
    apiGetUmap()
      .then((data) => {
        setPoints(data);
        setUsingMockData(false);
        setLoading(false);
      })
      .catch(() => {
        // If backend not available, use mock data
        setUsingMockData(true);
        setLoading(false);
      });
  }, []);

  return (
    <div style={{ fontFamily: "system-ui, Arial", background: "#fff", minHeight: "100vh" }}>
      <div style={{ width: "100%", maxWidth: 1200, margin: "0 auto", padding: "0 20px" }}>
        {/* TOP BAR */}
        <header
          style={{
            position: "sticky",
            top: 0,
            background: "white",
            borderBottom: "1px solid #e2e8f0",
            padding: "14px 0",
            zIndex: 10,
          }}
        >
          <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", gap: 16 }}>
            <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
              <span style={{ fontSize: 24 }}>🫀</span>
              <span style={{ fontWeight: 800, fontSize: 18 }}>hERG Cardiotoxicity Predictor</span>
            </div>

            <nav style={{ display: "flex", gap: 14, flexWrap: "wrap" }}>
              <a href="#about" style={navLink}>
                About
              </a>
              <a href="#demo" style={navLink}>
                Demo
              </a>
              <a href="#methodology" style={navLink}>
                Methodology
              </a>
              <a href="#dataset" style={navLink}>
                Dataset
              </a>
            </nav>

            <a href="#demo" style={primaryBtn}>
              Try Demo
            </a>
          </div>
        </header>

        {/* HERO */}
        <section style={{ padding: "48px 0 28px" }}>
          <h1 style={{ fontSize: 42, margin: 0, fontWeight: 900, lineHeight: 1.1 }}>
            Predict Cardiotoxicity of Molecules
          </h1>
          <p style={{ color: "#475569", marginTop: 12, fontSize: 17, maxWidth: 720, lineHeight: 1.6 }}>
            Analyze molecules for potential hERG channel blocking activity. Input a SMILES string or draw a molecule
            structure to get instant toxicity predictions with interpretable explanations.
          </p>

          <div style={{ display: "flex", gap: 10, marginTop: 16, flexWrap: "wrap" }}>
            <Badge text="Real-time Prediction" />
            <Badge text="Interpretable AI" />
            <Badge text="Chemical Space Visualization" />
            <Badge text="Lipinski Rule of 5" />
          </div>

          <div style={{ display: "flex", gap: 12, marginTop: 20, flexWrap: "wrap" }}>
            <a href="#demo" style={primaryBtn}>
              Start Analysis
            </a>
            <a href="#methodology" style={secondaryBtn}>
              Learn More
            </a>
          </div>
        </section>

        <Divider />

        {/* ABOUT */}
        <Section id="about" title="About">
          <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 24 }}>
            <div>
              <h3 style={{ margin: "0 0 8px", fontSize: 16, fontWeight: 700 }}>What is hERG Toxicity?</h3>
              <p style={{ margin: 0, lineHeight: 1.6, color: "#475569" }}>
                The hERG (human Ether-à-go-go-Related Gene) potassium channel plays a critical role in cardiac
                repolarization. Drug-induced hERG blockade can lead to QT interval prolongation and potentially fatal
                cardiac arrhythmias. Early prediction of hERG liability is essential in drug development.
              </p>
            </div>
            <div>
              <h3 style={{ margin: "0 0 8px", fontSize: 16, fontWeight: 700 }}>Our Approach</h3>
              <p style={{ margin: 0, lineHeight: 1.6, color: "#475569" }}>
                This tool uses machine learning to predict hERG cardiotoxicity from molecular structure. It combines
                MACCS fingerprints with physicochemical descriptors, trained on experimentally validated compounds. The
                model provides probability scores with interpretable explanations.
              </p>
            </div>
          </div>
        </Section>

        {/* DEMO - Main Feature */}
        <Section id="demo" title="Molecule Analysis Demo">
          {loading ? (
            <div style={{ padding: 40, textAlign: "center", color: "#64748b" }}>
              Loading chemical space data...
            </div>
          ) : (
            <>
              {usingMockData && (
                <div
                  style={{
                    padding: "8px 12px",
                    marginBottom: 12,
                    background: "#fef3c7",
                    border: "1px solid #f59e0b",
                    borderRadius: 8,
                    fontSize: 13,
                    color: "#92400e",
                  }}
                >
                  Using mock data - backend unavailable. Start the backend server to see real chemical space data.
                </div>
              )}
              <MoleculeAnalyzer backgroundPoints={points} />
            </>
          )}
        </Section>

        {/* METHODOLOGY */}
        <Section id="methodology" title="Methodology">
          <div style={{ display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 20 }}>
            <MethodCard
              icon="🧬"
              title="Feature Engineering"
              description="167 MACCS fingerprints capturing structural motifs combined with 130 physicochemical descriptors including LogP, molecular weight, and hydrogen bond counts."
            />
            <MethodCard
              icon="🤖"
              title="ML Model"
              description="XGBoost classifier trained on experimentally validated hERG binding data. Optimized for balanced sensitivity and specificity in toxicity prediction."
            />
            <MethodCard
              icon="📊"
              title="Visualization"
              description="UMAP dimensionality reduction projects molecules into 2D chemical space, allowing visualization of your molecule relative to known compounds."
            />
          </div>

          <div style={{ marginTop: 24 }}>
            <h3 style={{ margin: "0 0 12px", fontSize: 16, fontWeight: 700 }}>Explainability Features</h3>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(2, 1fr)", gap: 16 }}>
              <ExplainCard
                title="Structural Alerts"
                items={[
                  "Aromatic heterocycles (hERG binding risk)",
                  "Tertiary amines (ion channel interaction)",
                  "Quinone-like systems (oxidative stress)",
                ]}
              />
              <ExplainCard
                title="Physicochemical Heuristics"
                items={[
                  "High lipophilicity (LogP > 3.5)",
                  "High molecular weight (MW > 450 Da)",
                  "Multiple H-bond groups (HBD + HBA > 8)",
                ]}
              />
            </div>
          </div>
        </Section>

        {/* DATASET */}
        <Section id="dataset" title="Dataset">
          <div style={{ display: "grid", gridTemplateColumns: "repeat(4, 1fr)", gap: 16 }}>
            <StatCard value="~7,000" label="Compounds" />
            <StatCard value="Binary" label="Classification" />
            <StatCard value="297" label="Features" />
            <StatCard value="85%+" label="Accuracy" />
          </div>
          <p style={{ marginTop: 16, color: "#475569", lineHeight: 1.6 }}>
            The model is trained on publicly available hERG binding data from ChEMBL database, including both active
            blockers and inactive compounds. The dataset is balanced and curated to ensure reliable predictions.
          </p>
        </Section>

        <footer
          style={{
            padding: "32px 0",
            color: "#64748b",
            borderTop: "1px solid #e2e8f0",
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
          }}
        >
          <span>© {new Date().getFullYear()} hERG Cardiotoxicity Predictor</span>
          <span style={{ fontSize: 13 }}>Built with React + FastAPI + XGBoost</span>
        </footer>
      </div>
    </div>
  );
}

function Section(props: { id: string; title: string; children: ReactNode }) {
  return (
    <section id={props.id} style={{ padding: "34px 0" }}>
      <h2 style={{ margin: 0, fontSize: 26, fontWeight: 800 }}>{props.title}</h2>
      <div style={{ marginTop: 16, color: "#334155" }}>{props.children}</div>
      <Divider />
    </section>
  );
}

function Divider() {
  return <div style={{ borderTop: "1px solid #e2e8f0", marginTop: 24 }} />;
}

function Badge({ text }: { text: string }) {
  return (
    <span
      style={{
        padding: "6px 12px",
        borderRadius: 999,
        border: "1px solid #cbd5e1",
        background: "#f8fafc",
        fontSize: 13,
        fontWeight: 600,
        color: "#0f172a",
      }}
    >
      {text}
    </span>
  );
}

function MethodCard({ icon, title, description }: { icon: string; title: string; description: string }) {
  return (
    <div style={{ padding: 20, border: "1px solid #e2e8f0", borderRadius: 16, background: "#fafafa" }}>
      <div style={{ fontSize: 32, marginBottom: 12 }}>{icon}</div>
      <h3 style={{ margin: "0 0 8px", fontSize: 16, fontWeight: 700 }}>{title}</h3>
      <p style={{ margin: 0, fontSize: 14, color: "#64748b", lineHeight: 1.5 }}>{description}</p>
    </div>
  );
}

function ExplainCard({ title, items }: { title: string; items: string[] }) {
  return (
    <div style={{ padding: 16, border: "1px solid #e2e8f0", borderRadius: 12, background: "#fff" }}>
      <h4 style={{ margin: "0 0 10px", fontSize: 14, fontWeight: 700, color: "#0f172a" }}>{title}</h4>
      <ul style={{ margin: 0, paddingLeft: 20 }}>
        {items.map((item, i) => (
          <li key={i} style={{ fontSize: 13, color: "#64748b", lineHeight: 1.6 }}>
            {item}
          </li>
        ))}
      </ul>
    </div>
  );
}

function StatCard({ value, label }: { value: string; label: string }) {
  return (
    <div style={{ padding: 20, border: "1px solid #e2e8f0", borderRadius: 12, textAlign: "center" }}>
      <div style={{ fontSize: 28, fontWeight: 800, color: "#1d4ed8" }}>{value}</div>
      <div style={{ fontSize: 13, color: "#64748b", marginTop: 4 }}>{label}</div>
    </div>
  );
}

const navLink: CSSProperties = {
  color: "#475569",
  textDecoration: "none",
  fontWeight: 500,
  fontSize: 14,
};

const primaryBtn: CSSProperties = {
  background: "#1d4ed8",
  color: "white",
  padding: "10px 16px",
  borderRadius: 12,
  textDecoration: "none",
  fontWeight: 700,
  fontSize: 14,
};

const secondaryBtn: CSSProperties = {
  background: "white",
  color: "#0f172a",
  padding: "10px 16px",
  borderRadius: 12,
  textDecoration: "none",
  border: "1px solid #cbd5e1",
  fontWeight: 700,
  fontSize: 14,
};
