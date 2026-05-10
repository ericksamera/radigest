package gff

import "fmt"

const upperHex = "0123456789ABCDEF"

// EscapeSeqID percent-encodes a value for the first GFF3 column. This keeps
// ordinary accession/chromosome identifiers readable while escaping whitespace,
// separators, percent signs, and other characters that can break tab-delimited
// GFF3 parsing.
func EscapeSeqID(value string) string {
	if value == "" {
		return "."
	}
	return percentEncode(value, isSeqIDSafe)
}

// EscapeAttributeValue percent-encodes a GFF3 attribute value. In particular,
// it escapes ';', '=', ',', '%', tabs, and newlines so generated IDs and other
// values cannot create extra attributes or corrupt the ninth GFF3 column.
func EscapeAttributeValue(value string) string {
	return percentEncode(value, isAttributeValueSafe)
}

// FragmentAttributes builds the attributes used for radigest fragment features.
func FragmentAttributes(chr string, ordinal, length int) string {
	return fmt.Sprintf("ID=%s;Length=%d", fragmentID(chr, ordinal), length)
}

func fragmentID(chr string, ordinal int) string {
	if chr == "" {
		return fmt.Sprintf("frag%d", ordinal)
	}
	return fmt.Sprintf("%s_%d", EscapeAttributeValue(chr), ordinal)
}

func percentEncode(value string, safe func(byte) bool) string {
	needsEscape := false
	for i := 0; i < len(value); i++ {
		if !safe(value[i]) {
			needsEscape = true
			break
		}
	}
	if !needsEscape {
		return value
	}

	out := make([]byte, 0, len(value)+8)
	for i := 0; i < len(value); i++ {
		b := value[i]
		if safe(b) {
			out = append(out, b)
			continue
		}
		out = append(out, '%', upperHex[b>>4], upperHex[b&0x0f])
	}
	return string(out)
}

func isSeqIDSafe(b byte) bool {
	return isAlphaNum(b) || b == '.' || b == ':' || b == '^' || b == '*' || b == '$' || b == '@' || b == '!' || b == '+' || b == '_' || b == '?' || b == '-' || b == '|'
}

func isAttributeValueSafe(b byte) bool {
	return isAlphaNum(b) || b == '.' || b == ':' || b == '_' || b == '-' || b == '|'
}

func isAlphaNum(b byte) bool {
	return (b >= 'A' && b <= 'Z') || (b >= 'a' && b <= 'z') || (b >= '0' && b <= '9')
}
