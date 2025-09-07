#!/usr/bin/env python3
"""Example usage of the gene search functionality."""

import asyncio

from sirnaforge.data.gene_search import DatabaseType, GeneSearcher


async def main() -> None:
    """Example gene search demonstration."""

    # Initialize the gene searcher
    searcher = GeneSearcher(preferred_db=DatabaseType.ENSEMBL)

    # Example 1: Search for TP53 gene with metadata only (no sequences)
    print("🔍 Searching for TP53 gene (metadata only)...")
    result = await searcher.search_gene("TP53", include_sequence=False)

    if result.success and result.gene_info:
        print(f"✅ Found gene: {result.gene_info.gene_name} ({result.gene_info.gene_id})")
        print(f"📍 Location: {result.gene_info.chromosome}:{result.gene_info.start}-{result.gene_info.end}")
        print(f"🧬 Found {len(result.transcripts)} transcripts")

        # Show first few transcripts
        for i, transcript in enumerate(result.transcripts[:3]):
            canonical = "⭐" if transcript.is_canonical else ""
            print(f"   {i + 1}. {transcript.transcript_id} {canonical}")
    else:
        print(f"❌ Search failed: {result.error}")

    print("\n" + "=" * 60 + "\n")

    # Example 2: Search across multiple databases
    print("🔍 Searching BRCA1 across all databases...")
    results = await searcher.search_multiple_databases("BRCA1", include_sequence=False)

    for result in results:
        if result.success and result.gene_info:
            print(f"✅ {result.database.value}: {result.gene_info.gene_name} - {len(result.transcripts)} transcripts")
        else:
            print(f"❌ {result.database.value}: {result.error}")


if __name__ == "__main__":
    # Run the async example
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\n🛑 Search cancelled by user")
    except Exception as e:
        print(f"❌ Error: {e}")
