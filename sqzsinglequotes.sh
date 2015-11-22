# https://stackoverflow.com/questions/16881107/shell-script-rename-multi-files-and-remove-single-quote
find . -type f -name "*'*" | while IFS= read -r file
do
  # we need to avoid replacing characters in the path to the file,
  # so split it into dirname and filename.
  DIRNAME=$(dirname "$file")
  FILENAME=$(basename "$file")
  NEWNAME=$(sed "s/'//g" <<< "$FILENAME")
  mv -v --no-clobber "$file" "$DIRNAME/$NEWNAME" || echo "$DIRNAME/$NEWNAME already exists, not overwriting."
done
