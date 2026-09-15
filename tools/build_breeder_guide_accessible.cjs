#!/usr/bin/env node

"use strict";

const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");
const { chromium } = require("playwright");
const { PDFDocument } = require("pdf-lib");

const root = path.resolve(__dirname, "..");
const sourcePath = path.join(
  root,
  "inst",
  "guide",
  "HapBlockR_Breeder_Guide.md"
);
const outputPath = path.join(
  root,
  "inst",
  "extdata",
  "HapBlockR_Breeder_Guide.pdf"
);

function slugify(text) {
  return text
    .replace(/<[^>]+>/g, "")
    .replace(/&amp;/g, "and")
    .replace(/[^A-Za-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "")
    .toLowerCase();
}

function addHeadingIdentifiers(html) {
  const used = new Map();
  return html.replace(
    /<(h[1-3])>([\s\S]*?)<\/\1>/g,
    (match, tag, text) => {
      const base = slugify(text) || "section";
      const count = used.get(base) || 0;
      used.set(base, count + 1);
      const identifier = count ? `${base}-${count + 1}` : base;
      return `<${tag} id="${identifier}">${text}</${tag}>`;
    }
  );
}

function replaceContents(html) {
  const contents = `
<h2 id="contents">Contents</h2>
<nav class="contents" aria-label="Guide contents">
  <ol>
    <li><a href="#1-purpose-and-scope">Purpose and scope</a></li>
    <li><a href="#2-before-any-recommendation">Before any recommendation</a></li>
    <li><a href="#3-the-nine-decision-tools-and-their-variants">The nine decision tools and their variants</a></li>
    <li><a href="#4-worked-crossing-decision">Worked crossing decision</a></li>
    <li><a href="#5-interpreting-quality-control-and-uncertainty">Interpreting quality control and uncertainty</a></li>
    <li><a href="#6-interpreting-results-and-defining-their-scope">Interpreting results and defining their scope</a></li>
    <li><a href="#7-decision-sign-off">Decision sign-off</a></li>
    <li><a href="#8-supporting-tools-and-complete-function-map">Supporting tools and complete function map</a></li>
    <li><a href="#9-references">References and change history</a></li>
  </ol>
</nav>`;
  return html.replace(
    /<h2 id="contents">Contents<\/h2>[\s\S]*?(?=<h2 id="1-purpose-and-scope">)/,
    contents
  );
}

function documentHtml(body) {
  return `<!doctype html>
<html lang="en-GB">
<head>
  <meta charset="utf-8">
  <meta name="author" content="Félicien Akohoue">
  <meta name="description" content="HapBlockR breeder guide for traceable parent and cross decisions">
  <title>The HapBlockR Breeder's Guide</title>
  <style>
    :root {
      --navy: #123167;
      --blue: #2467a6;
      --green: #3d7d44;
      --ink: #172033;
      --muted: #536276;
      --line: #b8c8d8;
      --pale-blue: #edf5fb;
      --pale-green: #eef7ef;
    }
    @page {
      size: A4;
      margin: 17mm 16mm 20mm;
    }
    * { box-sizing: border-box; }
    html {
      font-family: "Arial", "Helvetica", sans-serif;
      color: var(--ink);
      font-size: 10.25pt;
      line-height: 1.38;
    }
    body { margin: 0; }
    h1, h2, h3 {
      color: var(--navy);
      line-height: 1.16;
      break-after: avoid;
    }
    h1 {
      margin: 0 0 5mm;
      padding-bottom: 4mm;
      border-bottom: 2.5pt solid var(--green);
      font-size: 25pt;
      letter-spacing: -0.35pt;
    }
    h2 {
      margin: 6mm 0 3mm;
      font-size: 16pt;
      border-bottom: 0.75pt solid var(--line);
      padding-bottom: 1.8mm;
    }
    h3 {
      margin: 4mm 0 2mm;
      font-size: 12.2pt;
      color: var(--green);
    }
    p { margin: 0 0 2.5mm; orphans: 3; widows: 3; }
    ul, ol { margin: 1.5mm 0 3mm 5mm; padding-left: 5mm; }
    li { margin: 0 0 1.1mm; }
    a { color: #0d5f91; text-decoration: underline; }
    code {
      font-family: "Consolas", "Courier New", monospace;
      font-size: 0.92em;
      background: #f2f4f7;
      border-radius: 2px;
      padding: 0.2mm 0.8mm;
    }
    pre {
      break-inside: avoid;
      white-space: pre-wrap;
      margin: 2mm 0 3mm;
    }
    pre code {
      display: block;
      padding: 2mm 2.5mm;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      margin: 2mm 0 4mm;
      font-size: 8.8pt;
      break-inside: avoid;
    }
    thead { display: table-header-group; }
    th {
      background: var(--navy);
      color: white;
      text-align: left;
      font-weight: 700;
      padding: 1.7mm 2mm;
      border: 0.5pt solid var(--navy);
    }
    th code {
      color: inherit;
      background: transparent;
      padding: 0;
    }
    td {
      vertical-align: top;
      padding: 1.5mm 2mm;
      border: 0.5pt solid var(--line);
    }
    tbody tr:nth-child(even) { background: var(--pale-blue); }
    blockquote {
      margin: 3mm 0;
      padding: 2.5mm 3.5mm;
      border-left: 3pt solid var(--green);
      background: var(--pale-green);
    }
    .page-break { break-before: page; height: 0; }
    .contents {
      background: var(--pale-blue);
      border: 0.75pt solid var(--line);
      border-radius: 4px;
      padding: 3mm 5mm;
    }
    .contents ol { margin-bottom: 0; }
    h1 + h2 {
      margin-top: 0;
      border: 0;
      color: var(--green);
      font-size: 15pt;
    }
    h1 + h2 + p,
    h1 + h2 + p + p,
    h1 + h2 + p + p + p {
      color: var(--muted);
    }
  </style>
</head>
<body>
<main aria-label="HapBlockR breeder guide">
${body}
</main>
</body>
</html>`;
}

(async () => {
  if (!fs.existsSync(sourcePath)) {
    throw new Error(`Guide source not found: ${sourcePath}`);
  }
  const { marked } = await import(
    pathToFileURL(require.resolve("marked")).href
  );
  let markdown = fs.readFileSync(sourcePath, "utf8");
  markdown = markdown.replace(
    /<!--\s*pagebreak\s*-->/gi,
    '<div class="page-break" aria-hidden="true"></div>'
  );
  let body = marked.parse(markdown, { gfm: true });
  body = addHeadingIdentifiers(body);
  body = replaceContents(body);

  const browserCandidates = [
    process.env.CHROME_PATH,
    chromium.executablePath(),
    "C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe",
    "C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe",
    "/usr/bin/google-chrome",
    "/usr/bin/chromium",
    "/usr/bin/chromium-browser"
  ].filter(Boolean);
  const browserExecutable = browserCandidates.find(fs.existsSync);
  if (!browserExecutable) {
    throw new Error(
      "No Chromium executable was found. Set CHROME_PATH or install the " +
      "Playwright Chromium browser."
    );
  }
  const browser = await chromium.launch({
    headless: true,
    executablePath: browserExecutable
  });
  try {
    const page = await browser.newPage();
    await page.setContent(documentHtml(body), { waitUntil: "load" });
    await page.pdf({
      path: outputPath,
      format: "A4",
      printBackground: true,
      preferCSSPageSize: true,
      tagged: true,
      outline: true,
      displayHeaderFooter: true,
      headerTemplate: "<span></span>",
      footerTemplate: `
        <div style="width:100%;font:8px Arial;color:#536276;
                    padding:0 16mm;display:flex;justify-content:space-between">
          <span>HBR-GUIDE-001 · HapBlockR 0.3.12.9000</span>
          <span>Page <span class="pageNumber"></span> of
                <span class="totalPages"></span></span>
        </div>`
    });
  } finally {
    await browser.close();
  }

  const pdf = await PDFDocument.load(fs.readFileSync(outputPath));
  pdf.setTitle("The HapBlockR Breeder's Guide");
  pdf.setAuthor("Félicien Akohoue");
  pdf.setSubject("Traceable parent and cross decisions with HapBlockR");
  pdf.setKeywords([
    "HapBlockR",
    "breeding",
    "parent selection",
    "genomic mating"
  ]);
  pdf.setCreator("HapBlockR accessible guide builder");
  fs.writeFileSync(outputPath, await pdf.save());

  const bytes = fs.statSync(outputPath).size;
  if (bytes < 20000) {
    throw new Error(`Generated guide is unexpectedly small: ${bytes} bytes`);
  }
  process.stdout.write(`Generated ${outputPath} (${bytes} bytes)\n`);
})().catch((error) => {
  process.stderr.write(`${error.stack || error.message}\n`);
  process.exit(1);
});
