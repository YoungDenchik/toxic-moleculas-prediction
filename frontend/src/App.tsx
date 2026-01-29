import { useEffect, useState } from "react";
import type { ReactNode } from "react";
import type { CSSProperties } from "react";


import { apiGetUmap } from "./api/client";
import { UmapPlot } from "./components/UmapPlot";
import { PredictPanel } from "./components/PredictPanel.tsx";
import type { UmapPoint } from "./types";
import { mockUmap } from "./mocks/umap";

export default function App() {
  const [points, setPoints] = useState<UmapPoint[]>(mockUmap);
  const [picked, setPicked] = useState<UmapPoint | null>(null);

  useEffect(() => {
    // коли бек буде готовий — підтягуємо з API
    apiGetUmap()
      .then(setPoints)
      .catch(() => {
        // якщо бек не працює — лишаємось на mock
      });
  }, []);

  return (
    <div style={{ fontFamily: "system-ui, Arial", background: "#fff", minHeight: "100vh" }}>
      <div style={{ width: "100%", maxWidth: 1100, margin: "0 auto", padding: "0 16px" }}>
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
            <div style={{ fontWeight: 800, fontSize: 18 }}>hERG Predictor</div>

            <nav style={{ display: "flex", gap: 14, flexWrap: "wrap" }}>
              <a href="#about">About</a>
              <a href="#demo">Demo</a>
              <a href="#tracks">Tracks</a>
              <a href="#dataset">Dataset</a>
              <a href="#results">Results</a>
            </nav>

            <a
              href="#demo"
              style={{
                background: "#1d4ed8",
                color: "white",
                padding: "10px 12px",
                borderRadius: 12,
                textDecoration: "none",
                fontWeight: 600,
              }}
            >
              Try Demo
            </a>
          </div>
        </header>

        {/* HERO */}
        <section style={{ padding: "48px 0 28px" }}>
          <h1 style={{ fontSize: 42, margin: 0, fontWeight: 900, lineHeight: 1.1 }}>
            Predict hERG toxicity from molecules
          </h1>
          <p style={{ color: "#475569", marginTop: 12, fontSize: 16, maxWidth: 720 }}>
            Web app that takes a molecule (SMILES) and returns hERG risk + explanation.
          </p>

          <div style={{ display: "flex", gap: 10, marginTop: 14, flexWrap: "wrap" }}>
            <Badge text="Real-time" />
            <Badge text="Explanation" />
          </div>

          <div style={{ display: "flex", gap: 12, marginTop: 18, flexWrap: "wrap" }}>
            <a href="#demo" style={primaryBtn}>
              Start demo
            </a>
            <a href="#tracks" style={secondaryBtn}>
              See tracks
            </a>
          </div>
        </section>

        <Divider />

        {/* SECTIONS */}
        <Section id="about" title="About">
          Тут буде короткий опис проєкту (1–2 абзаци).
        </Section>

        <Section id="demo" title="Demo">
          <div style={{ marginTop: 14 }}>
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16 }}>
              <PredictPanel />
              <UmapPlot points={points} onPick={setPicked} />
            </div>

            {picked && (
              <div style={{ marginTop: 12, color: "#334155" }}>
                <b>Selected:</b> {picked.id} {picked.smiles ? `• ${picked.smiles}` : ""}
              </div>
            )}
          </div>
        </Section>

        <Section id="tracks" title="Tracks">
          Тут буде Track A/B/C (коротко).
        </Section>

        <Section id="dataset" title="Dataset">
          Тут буде опис датасету: розмір, що за мітки, звідки.
        </Section>

        <Section id="results" title="Results">
          Тут будуть метрики, таблиця/графіки.
        </Section>

        <Section id="mentor" title="Mentor">
          Тут буде ментор/команда.
        </Section>

        <Section id="links" title="Links">
          Тут будуть посилання на GitHub/презентацію/демо.
        </Section>

        <footer style={{ padding: "32px 0", color: "#64748b", borderTop: "1px solid #e2e8f0" }}>
          © {new Date().getFullYear()} ML Week Project
        </footer>
      </div>
    </div>
  );
}

function Section(props: { id: string; title: string; children: ReactNode }) {
  return (
    <section id={props.id} style={{ padding: "34px 0" }}>
      <h2 style={{ margin: 0, fontSize: 26, fontWeight: 800 }}>{props.title}</h2>
      <div style={{ marginTop: 10, color: "#334155" }}>{props.children}</div>
      <Divider />
    </section>
  );
}

function Divider() {
  return <div style={{ borderTop: "1px solid #e2e8f0" }} />;
}

function Badge({ text }: { text: string }) {
  return (
    <span
      style={{
        padding: "6px 10px",
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

const primaryBtn: CSSProperties = {
  background: "#1d4ed8",
  color: "white",
  padding: "10px 14px",
  borderRadius: 12,
  textDecoration: "none",
  fontWeight: 700,
};

const secondaryBtn: CSSProperties = {
  background: "white",
  color: "#0f172a",
  padding: "10px 14px",
  borderRadius: 12,
  textDecoration: "none",
  border: "1px solid #cbd5e1",
  fontWeight: 700,
};
