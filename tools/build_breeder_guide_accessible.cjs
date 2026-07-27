#!/usr/bin/env node

"use strict";

const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");
const { chromium } = require("playwright");
const { PDFDocument } = require("pdf-lib");

const root = path.resolve(__dirname, "..");
// The version-controlled source is a genuine R Markdown document (YAML
// frontmatter, one inline `r ...` expression for the build date). This
// script does not invoke R/knitr/pandoc -- it parses the document as plain
// Markdown via `marked` -- so parseFrontmatter()/stripHeadingAttributes()/
// numberHeadings() below reproduce, by hand, the small subset of
// rmarkdown::html_document behaviour (title block, `number_sections`,
// `{.unnumbered}`) that this guide actually relies on.
const sourcePath = path.join(
  root,
  "inst",
  "guide",
  "HapBlockR_Breeder_Guide.Rmd"
);
const pdfOutputPath = path.join(
  root,
  "inst",
  "extdata",
  "HapBlockR_Breeder_Guide.pdf"
);
const htmlOutputPath = path.join(
  root,
  "inst",
  "extdata",
  "HapBlockR_Breeder_Guide.html"
);

function slugify(text) {
  return text
    .replace(/<[^>]+>/g, "")
    .replace(/&amp;/g, "and")
    .replace(/[^A-Za-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "")
    .toLowerCase();
}

// Extracts the leading `---\n...\n---` YAML block by hand (no YAML
// dependency is otherwise needed), evaluates the one inline R date
// expression at build time, and returns the frontmatter fields plus the
// markdown body with the block removed.
function parseFrontmatter(markdown) {
  const match = markdown.match(/^---\n([\s\S]*?)\n---\n?/);
  if (!match) {
    throw new Error("Expected YAML frontmatter at the top of the guide source.");
  }
  const block = match[1];
  const field = (name) => {
    const m = block.match(new RegExp(`^${name}:\\s*"(.*)"\\s*$`, "m"));
    return m ? m[1] : "";
  };
  let date = field("date");
  if (/`r .*`/.test(date)) {
    date = new Date().toLocaleDateString("en-GB", {
      day: "numeric",
      month: "long",
      year: "numeric"
    });
  }
  const meta = {
    title: field("title"),
    subtitle: field("subtitle"),
    author: field("author"),
    date,
    compatibleVersion: field("compatible_package_version")
  };
  return { meta, body: markdown.slice(match[0].length) };
}

// Removes pandoc header-attribute syntax (`{.unnumbered}`, `{#id}`, ...)
// from ATX headings before handing the text to `marked`, which does not
// understand it. `.unnumbered` is preserved as an inline HTML comment so
// numberHeadings() can still see it after parsing.
function stripHeadingAttributes(markdown) {
  return markdown.replace(
    /^(#{1,6}\s+.*?)\s*\{([^}]*)\}\s*$/gm,
    (match, heading, attrs) => (
      /\.unnumbered/.test(attrs) ? `${heading} <!--unnumbered-->` : heading
    )
  );
}

// Reproduces pandoc's `number_sections: true`: walks h1/h2/h3 in document
// order, numbers them "1", "1.1", "1.1.1" (skipping any marked
// `<!--unnumbered-->`), and assigns each a slug id built from the now
// numbered text -- so `1 Purpose and scope` gets id="1-purpose-and-scope",
// matching how a real pandoc/knitr render would build the same guide.
// Returns the updated HTML plus the ordered list of numbered h1 entries,
// which insertContents() below uses to build a Contents block that can
// never drift out of sync with the actual chapters.
function numberHeadings(html) {
  const counters = [0, 0, 0];
  const used = new Map();
  const h1Entries = [];
  const updated = html.replace(
    /<(h[1-3])>([\s\S]*?)<\/\1>/g,
    (match, tag, rawText) => {
      const level = Number(tag[1]);
      const unnumbered = /<!--unnumbered-->/.test(rawText);
      const text = rawText.replace(/\s*<!--unnumbered-->\s*/g, "").trim();
      let displayText = text;
      if (!unnumbered) {
        counters[level - 1] += 1;
        for (let i = level; i < 3; i++) counters[i] = 0;
        const prefix = counters.slice(0, level).join(".");
        displayText = `${prefix} ${text}`;
      }
      const base = slugify(displayText) || "section";
      const count = used.get(base) || 0;
      used.set(base, count + 1);
      const id = count ? `${base}-${count + 1}` : base;
      if (level === 1 && !unnumbered) {
        h1Entries.push({ id, text });
      }
      return `<${tag} id="${id}">${displayText}</${tag}>`;
    }
  );
  return { html: updated, h1Entries };
}

// Inserts a Contents block built from the real numbered chapters, right
// before the first one, instead of relying on a hand-maintained anchor
// list that silently drifts out of date whenever a chapter is added,
// removed or renamed.
function insertContents(html, h1Entries) {
  if (h1Entries.length === 0) return html;
  const items = h1Entries
    .map((entry) => `    <li><a href="#${entry.id}">${entry.text}</a></li>`)
    .join("\n");
  const contents = `
<h2 id="contents">Contents</h2>
<nav class="contents" aria-label="Guide contents">
  <ol>
${items}
  </ol>
</nav>
`;
  const anchor = `<h1 id="${h1Entries[0].id}">`;
  const idx = html.indexOf(anchor);
  if (idx === -1) return html;
  return html.slice(0, idx) + contents + html.slice(idx);
}

// Builds the pandoc-style title block (title/subtitle/author/date) that
// rmarkdown::html_document would normally generate from the YAML
// frontmatter, since `marked` never sees the frontmatter at all. The CSS in
// documentHtml() below (`h1 + h2`, `h1 + h2 + p ...`) targets exactly this
// four-element structure.
function titleBlockHtml(meta) {
  const compat = meta.compatibleVersion
    ? `<p class="compat">Compatible with HapBlockR ${meta.compatibleVersion}</p>`
    : "";
  return `<h1 class="title">${meta.title}</h1>
<h2 class="subtitle">${meta.subtitle}</h2>
<p class="author">${meta.author}</p>
<p class="date">${meta.date}</p>
${compat}`;
}

function documentHtml(body, meta) {
  return `<!doctype html>
<html lang="en-GB">
<head>
  <meta charset="utf-8">
  <meta name="author" content="${meta.author}">
  <meta name="description" content="HapBlockR breeder guide for traceable parent and cross decisions">
  <title>${meta.title}</title>
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
    /* Title-block elements are matched by class, not by h1+h2 tag
       adjacency: several real chapters (e.g. "Before any recommendation")
       have their first h2 subsection immediately after the chapter h1 with
       no intervening paragraph, so a bare adjacent-sibling selector would
       misapply the subtitle style to those subsection headings too. */
    h2.subtitle {
      margin-top: 0;
      border: 0;
      color: var(--green);
      font-size: 15pt;
    }
    p.author, p.date, p.compat {
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
  const raw = fs.readFileSync(sourcePath, "utf8");
  const { meta, body: withoutFrontmatter } = parseFrontmatter(raw);
  let markdown = stripHeadingAttributes(withoutFrontmatter);
  markdown = markdown.replace(
    /<!--\s*pagebreak\s*-->/gi,
    '<div class="page-break" aria-hidden="true"></div>'
  );
  let body = marked.parse(markdown, { gfm: true });
  const numbered = numberHeadings(body);
  body = insertContents(numbered.html, numbered.h1Entries);
  body = `${titleBlockHtml(meta)}\n${body}`;

  const fullHtml = documentHtml(body, meta);
  fs.writeFileSync(htmlOutputPath, fullHtml);
  process.stdout.write(`Generated ${htmlOutputPath}\n`);

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
    await page.setContent(fullHtml, { waitUntil: "load" });
    await page.pdf({
      path: pdfOutputPath,
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
          <span>HBR-GUIDE-001 · HapBlockR ${meta.compatibleVersion}</span>
          <span>Page <span class="pageNumber"></span> of
                <span class="totalPages"></span></span>
        </div>`
    });
  } finally {
    await browser.close();
  }

  const pdf = await PDFDocument.load(fs.readFileSync(pdfOutputPath));
  pdf.setTitle(meta.title);
  pdf.setAuthor(meta.author);
  pdf.setSubject("Traceable parent and cross decisions with HapBlockR");
  pdf.setKeywords([
    "HapBlockR",
    "breeding",
    "parent selection",
    "genomic mating"
  ]);
  pdf.setCreator("HapBlockR accessible guide builder");
  fs.writeFileSync(pdfOutputPath, await pdf.save());

  const bytes = fs.statSync(pdfOutputPath).size;
  if (bytes < 20000) {
    throw new Error(`Generated guide is unexpectedly small: ${bytes} bytes`);
  }
  process.stdout.write(`Generated ${pdfOutputPath} (${bytes} bytes)\n`);
})().catch((error) => {
  process.stderr.write(`${error.stack || error.message}\n`);
  process.exit(1);
});
