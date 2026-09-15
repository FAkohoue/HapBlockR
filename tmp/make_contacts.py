from pathlib import Path
import re

from PIL import Image, ImageDraw


folder = Path("tmp/guide-docx-pages-verified")
output = folder / "contacts"
output.mkdir(exist_ok=True)
pages = sorted(
    folder.glob("page-*.png"),
    key=lambda path: int(re.search(r"([0-9]+)$", path.stem).group(1)),
)

for start in range(0, len(pages), 4):
    batch = pages[start : start + 4]
    images = [Image.open(path).convert("RGB") for path in batch]
    width = max(image.width for image in images)
    height = max(image.height for image in images)
    sheet = Image.new("RGB", (2 * width, 2 * height + 40), "white")
    for index, image in enumerate(images):
        sheet.paste(image, ((index % 2) * width, (index // 2) * height + 40))
    page_numbers = [
        str(int(re.search(r"([0-9]+)$", path.stem).group(1))) for path in batch
    ]
    ImageDraw.Draw(sheet).text(
        (10, 10), "Pages " + ", ".join(page_numbers), fill="black"
    )
    sheet.save(output / f"contact-{start // 4 + 1:02d}.png")
