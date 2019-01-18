[![Build Status](https://travis-ci.org/mklarqvist/tomahawk.svg?branch=master)](https://travis-ci.org/mklarqvist/tomahawk)
[![Release](https://img.shields.io/badge/Release-beta_0.7.0-blue.svg)](https://github.com/mklarqvist/tomahawk/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<div align="center">
<img src="images/tomahawk.png" style="max-width:400px;">
</div>

# Fast calculation of LD in large-scale cohorts
Tomahawk is a machine-optimized library for computing
[linkage-disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium)
from population-sized datasets. Tomahawk permits close to real-time analysis of
regions-of-interest in datasets of many millions of diploid individuals on a
standard laptop. All algorithms are embarrassingly parallel and have been
successfully tested on chromosome-sized datasets with up to _10 million_
individuals.

Tomahawk uniquely constructs complete haplotype/genotype contigency matrices for
each comparison, perform statistical tests on the output data, and provide a
framework for querying the resulting data.

## CLI Commands

| Command        | Description                                                 |
|----------------|-------------------------------------------------------------|
| [`aggregate`](cli/cli-aggregate)    | data rasterization framework for `TWO` files                |
| `calc`         | calculate linkage disequilibrium                            |
| `concat`       | concatenate `TWO` files from the same set of samples        |
| `import`       | import `VCF`/`VCF.gz`/`BCF` to `TWK`                        |
| `sort`         | sort `TWO` file                                             |
| `view`         | `TWO`-&gt;`LD`/`TWO` view, `TWO` subset and filter          |
| `haplotype`    | extract per-sample haplotype strings in `FASTA`/binary format |
| `relationship` | compute marker-based pair-wise sample relationship matrices |
| `decay`        | compute LD-decay over distance                              |
| `prune`        | perform graph-based LD-pruning of variant sites             |

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.

#### Localization

> Default: `en`

Material for MkDocs supports internationalization (i18n) and provides
translations for all template variables and labels in the following languages:

<table style="white-space: nowrap;">
  <thead>
    <tr>
      <th colspan="4">Available languages</td>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>ar</code> / Arabic</td>
      <td><code>ca</code> / Catalan</td>
      <td><code>cs</code> / Czech</td>
      <td><code>da</code> / Danish</td>
    </tr>
    <tr>
      <td><code>nl</code> / Dutch</td>
      <td><code>en</code> / English</td>
      <td><code>fi</code> / Finnish</td>
      <td><code>fr</code> / French</td>
    </tr>
    <tr>
      <td><code>gl</code> / Galician</td>
      <td><code>de</code> / German</td>
      <td><code>he</code> / Hebrew</td>
      <td><code>hi</code> / Hindi</td>
    </tr>
    <tr>
      <td><code>hr</code> / Croatian</td>
      <td><code>hu</code> / Hungarian</td>
      <td><code>id</code> / Indonesian</td>
      <td><code>it</code> / Italian</td>
    </tr>
    <tr>
      <td><code>ja</code> / Japanese</td>
      <td><code>kr</code> / Korean</td>
      <td><code>no</code> / Norwegian</td>
      <td><code>fa</code> / Persian</td>
    </tr>
    <tr>
      <td><code>pl</code> / Polish</td>
      <td><code>pt</code> / Portugese</td>
      <td><code>ru</code> / Russian</td>
      <td><code>sr</code> / Serbian</td>
    </tr>
    <tr>
      <td><code>sh</code> / Serbo-Croatian</td>
      <td><code>sk</code> / Slovak</td>
      <td><code>es</code> / Spanish</td>
      <td><code>sv</code> / Swedish</td>
    </tr>
    <tr>
      <td colspan="2">
        <code>zh</code> / Chinese (Simplified)
      </td>
      <td colspan="2">
        <code>zh-Hant</code> / Chinese (Traditional)
      </td>
    </tr>
    <tr>
      <td colspan="2">
        <code>zh-TW</code> / Chinese (Taiwanese)
      </td>
      <td><code>tr</code> / Turkish</td>
      <td><code>uk</code> / Ukrainian</td>
    </tr>
    <tr>
      <td><code>vi</code> / Vietnamese</td>
      <td colspan="3" align="right">
        <a href="http://bit.ly/2EbzFc8">Submit a new language</a>
      </td>
    </tr>
  </tbody>
</table>

Specify the language with:

``` yaml
theme:
  language: 'en'
```

If the language is not specified, Material falls back to English. To create a
translation for another language, copy the localization file of an existing
language, name the new file using the [2-letter language code][16] and adjust
all translations:

``` sh
cp partials/language/en.html partials/language/jp.html
```

  [16]: https://www.w3schools.com/tags/ref_language_codes.asp

!!! success "Installation complete"
    
    Success: Tomahawk is now installed. Read the [tutorials](tutorial) to get started.

!!! warning "MkDocs 1.0 compatibility"

    While MkDocs 1.0 supports prebuilding the search index, Material currently
    doesn't support this setting as the default search behavior of the original
    theme was heavily modified for the sake of a better UX. Integration is
    possible, but a small subset of the features Material provides will not be
    portable to the prebuilt index mainly due to missing localization.

!!! warning "Only specify the languages you really need"

    Be aware that including support for other languages increases the general
    JavaScript payload by around 20kb (without gzip) and by another 15-30kb per
    language.

The separator for tokenization can be customized which makes it possible

!!! failure "Error: unrecognized theme 'material'"

    If you run into this error, the most common reason is that you installed
    MkDocs through some package manager (e.g. Homebrew or `apt-get`) and the
    Material theme through `pip`, so both packages end up in different
    locations. MkDocs only checks its install location for themes.

``` python hl_lines="3 4"
""" Bubble sort """
def bubble_sort(items):
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
```


<div align="center">
<img src="images/tomahawk_overview_problem.jpg">
</div>

<div align="center">
<img src="images/ld_overview.jpg">
</div>