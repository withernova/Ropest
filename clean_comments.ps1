# 清理 C# 文件中的大段注释
# 规则：
# 1. 保留 TODO/FIXME/HACK/NOTE/BUG 等功能性标记注释
# 2. 保留 XML 文档注释 (///)
# 3. 删除大段多行注释 (/* ... */)
# 4. 删除连续3行以上的单行注释块
# 5. 删除 Unity 自动生成的模板注释
# 6. 删除行尾注释中解释显而易见的内容（但保留简短必要的）

$scriptsPath = "e:\tmp_file_from_G\ROD\Assets\Scripts"
$csFiles = Get-ChildItem -Path $scriptsPath -Recurse -Filter "*.cs"

$unityTemplateComments = @(
    "Start is called once before the first execution of Update after the MonoBehaviour is created",
    "Update is called once per frame",
    "Start is called before the first frame update",
    "Called once per frame"
)

$functionalMarkers = @("TODO", "FIXME", "HACK", "NOTE", "BUG", "WARNING", "IMPORTANT")

foreach ($file in $csFiles) {
    $content = Get-Content -Path $file.FullName -Raw -Encoding UTF8
    $lines = Get-Content -Path $file.FullName -Encoding UTF8
    $newLines = @()
    $skipCount = 0
    
    for ($i = 0; $i -lt $lines.Count; $i++) {
        if ($skipCount -gt 0) {
            $skipCount--
            continue
        }
        
        $line = $lines[$i]
        $trimmed = $line.Trim()
        
        # 保留 XML 文档注释
        if ($trimmed.StartsWith("///")) {
            $newLines += $line
            continue
        }
        
        # 检查是否是功能性标记注释 (TODO/FIXME等)
        $isFunctional = $false
        foreach ($marker in $functionalMarkers) {
            if ($trimmed -match "//\s*$marker" -or $trimmed -match "/\*\s*$marker") {
                $isFunctional = $true
                break
            }
        }
        if ($isFunctional) {
            $newLines += $line
            continue
        }
        
        # 删除多行注释 /* ... */
        if ($trimmed.StartsWith("/*")) {
            # 检查是否在单行内结束
            if ($trimmed.Contains("*/")) {
                # 单行多行注释，删除
                continue
            }
            # 跳过多行注释
            $skipCount = 0
            for ($j = $i + 1; $j -lt $lines.Count; $j++) {
                $skipCount++
                if ($lines[$j].Contains("*/")) {
                    break
                }
            }
            continue
        }
        
        # 删除 Unity 模板注释
        $isTemplate = $false
        foreach ($template in $unityTemplateComments) {
            if ($trimmed -match [regex]::Escape($template)) {
                $isTemplate = $true
                break
            }
        }
        if ($isTemplate) {
            continue
        }
        
        # 检查连续注释块
        if ($trimmed.StartsWith("//") -and -not $trimmed.StartsWith("///")) {
            # 向前看，统计连续注释行数
            $consecutiveComments = 1
            for ($j = $i + 1; $j -lt $lines.Count; $j++) {
                $nextTrimmed = $lines[$j].Trim()
                if ($nextTrimmed.StartsWith("//") -and -not $nextTrimmed.StartsWith("///")) {
                    # 检查是否是功能性注释
                    $isFunc = $false
                    foreach ($marker in $functionalMarkers) {
                        if ($nextTrimmed -match "//\s*$marker") {
                            $isFunc = $true
                            break
                        }
                    }
                    if ($isFunc) { break }
                    $consecutiveComments++
                } elseif ([string]::IsNullOrWhiteSpace($nextTrimmed)) {
                    continue
                } else {
                    break
                }
            }
            
            # 如果是3行以上的注释块，跳过全部
            if ($consecutiveComments -ge 3) {
                $skipCount = $consecutiveComments - 1
                continue
            }
        }
        
        # 删除行尾的显而易见注释（但保留在代码后的简短注释）
        # 如果一行中有代码 + // 注释，且注释内容过于冗长或显而易见，删除注释部分
        if ($line -match "^(.*?)//(.+)$" -and -not $line.Trim().StartsWith("//")) {
            $codePart = $matches[1]
            $commentPart = $matches[2].Trim()
            
            # 保留短注释（少于15个字或包含关键信息）
            if ($commentPart.Length -gt 15) {
                $newLines += $codePart.TrimEnd()
                continue
            }
        }
        
        $newLines += $line
    }
    
    # 写回文件
    $newContent = $newLines -join "`n"
    # 确保文件以换行符结尾
    if (-not $newContent.EndsWith("`n")) {
        $newContent += "`n"
    }
    Set-Content -Path $file.FullName -Value $newContent -Encoding UTF8 -NoNewline
    Write-Host "Cleaned: $($file.FullName)"
}

Write-Host "`nDone! Cleaned $($csFiles.Count) files."
