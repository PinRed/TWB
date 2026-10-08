git remote add origin https://github.com/PinRed/TWB

find . -name "*.sh" -o -name "*.R" | xargs git add

git commit -m "Add all .sh and .R file"

git push -u origin main 
